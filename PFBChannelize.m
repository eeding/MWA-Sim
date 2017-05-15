function output = PFBChannelize(FS,S_IN,CONFIG,CHSEL,CHGAIN,QB,BITHIST,DOFLOAT)
%PFBCHANNELIZE uses PFBFilterNoSynth and radix2fft to perform the
%              two stage PFB filtering in the MWA DSP simulation.
%   PFBChannelize(FS,S_IN,CONFIG,CHSEL,CHGAIN,QB,BITHIST) - passes S_IN
%   through the 2-stage PFB simulation and produces an output containing 
%   the signal at various stages of the process.
%   
%   -> S_IN must be a single row vector containing the input signal. S_IN
%      can be floating point, but it will be implicitly casted into fi
%      representation in the filtering process.
%   -> CONFIG must be a list of two structs: 
%      {1} for the first PFB, {2} for the second PFB
%      Each struct must contain the following fields:
%      coeff - The floating point coefficients for the PFB
%      fi_coeff - The fixed point coefficients for the PFB
%      twiddle - The fixed point twiddle factors for the FFT in the PFB
%      output_nt - The Numeric Type requirement for the output of the PFB.
%                  (See MATLAB Documentation for information regarding FI
%                  objects and Numeric Types).
%   -> CHSEL must be a row vector containing the channel indices of the
%      coarse channels selected to continue onto the second PFB. Index 1 is
%      baseband, Index 256 is highest subband. Invalid values will
%      automatically be corrected in a rounding modulo 256 scheme, but
%      behavior is not guaranteed.
%   -> QB must be a vector of length 2 containing information about the
%      number of bits kept in the re-quantization that occcurs after each
%      stage. To match the physical system, use [5,4], which makes 5 bit
%      quantization after PFB1 and 4 bit quantization after PFB2.
%      Constraints: 0 < QB(1) < 15, 0 < QB(2) < 23, integer values only.
%      Set any (or both) QB value to 0 to skip re-quantization step for the
%      corresponding stage.
%   -> BITHIST is a debug flag. If set, bit histograms at various points
%      will be computed and included in the output.
%   -> DOFLOAT is a debug flag. If set, all computation will be performed
%      using floating point.
%
%   -> Outputs a list containing two structs:
%      {1} has data for first stage, {2} has data for the second stage
%      Each struct should contain the following fields:
%      out - The output of the PFB after the re-quantization
%      fbins - The frequency values associated with the channels in the
%              output.
%      (More outputs TBD)
%      If BITHIST is set, there will be additional fields with name prefix
%      "bithist" that correspond to the bit histograms at various points in
%      the particular PFB stage.


    output = {struct(),struct()};
    % Set maximum number of bits for bit histogram if applicable
    if ~DOFLOAT && BITHIST
        bmax = 32;
    end

    n1 = size(CONFIG{1}.fi_coeff,1); % this should be 512
    n2 = size(CONFIG{2}.coeff,1); % this should be 128
    if (~isequal(size(size(CHSEL)),[1,2]))
        error('ERROR: CHSEL must be row vector!');
    else
        CHSEL = mod(abs(CHSEL-1),ceil(n1/2))+1;
    end
    if ~isequal(size(QB),[1,2])
        error('ERROR: Need 1 x 2 matrix indicating quantization bit levels.');
    end
    if ~isequal(size(CHGAIN),size(CHSEL))
        disp('WARNING: CHGAIN dimension does not match CHSEL dimension. Assuming unity gain.');
        CHGAIN = ones(size(CHSEL));
    end
   
    % Perform first stage filtering
    if DOFLOAT
        X = PFBfiltNoSynth(double(S_IN),CONFIG{1}.coeff,0);
    else
        X = PFBfiltNoSynth(S_IN,CONFIG{1}.fi_coeff,CONFIG{1}.output_nt);
    end
    if size(X,2) < 1
        error('Insufficient data for 1st FFT.');
    end
    if DOFLOAT
        spec = radix2fft(X,radix2twiddles(n1));
    else
        spec = radix2fft(X,CONFIG{1}.twiddle);
    end
    freq = [(0:1/n1:0.5),((-0.5+1/n1):1/n1:-1/n1)].';
    freq = FS.*freq;

    % Quantize to fewer bits if fixed point and QB(1) is set
    if ~DOFLOAT && (QB(1) > 0 && QB(1) < CONFIG{1}.output_nt.WordLength-1)
        spec = quantize(spec/2^(16-QB(1)),numerictype(1,QB(1),0),'Round','Saturate');
    end
    
    % Store results for first stage
    if DOFLOAT
        output{1}.out = spec;
    else
        output{1}.out = spec.removefimath;
    end
    % Only take first half of spectrum for PFB1
    output{1}.out = output{1}.out(1:ceil(n1/2),:,:);
    output{1}.fbins = freq(1:ceil(n1/2),:,:);
    
    % Calculate the bit histograms for stage 1
    if exist('bmax','var')
        output{1}.bithist_inp = fiBitHist(S_IN,bmax);
        output{1}.bithist_filtout = fiBitHist(X,bmax);
        output{1}.bithist_fftout = fiBitHist(spec,bmax);
    end
    
    % Check if there's enough datapoints to continue to stage 2
    if(size(output{1}.out,2) < length(CONFIG{2}.coeff(:)))
        % Returning instead of throwing error so first stage data is still
        % returned
        disp("Not enough samples to proceed to stage 2!");
        return
    else
        % Select only a few coarse channels for next stage
        S_IN = output{1}.out(CHSEL,:,:);
        coarse_freq = output{1}.fbins(CHSEL);
        % Scale each channel by appropriate gain
        S_IN = S_IN.*repmat(CHGAIN.',1,size(S_IN,2),size(S_IN,3));
        % Move coarse channel index to dimension 3 so they will be acted on
        % in parallel for second stage
        X = permute(S_IN,[3,2,1]);
        
        % Perform second stage filtering
        if DOFLOAT
            X = PFBfiltNoSynth(X,CONFIG{2}.coeff,0);
            spec = radix2fft(X,radix2twiddles(n2));
        else
            X = PFBfiltNoSynth(X,CONFIG{2}.fi_coeff,CONFIG{2}.output_nt);
            spec = radix2fft(X,CONFIG{2}.twiddle);
        end
        % 2nd stage frequency bin is relative to each coarse bin
        freq = [(0:1/n2:0.5),((-0.5+1/n2):1/n2:-1/n2)].';
        freq = (FS/n1).*freq;
        % Reorder to make sure frequencies come out negative first
        spec = [spec(freq<0,:,:);spec(freq>=0,:,:)];
        freq = [freq(freq<0);freq(freq>=0)];
        
        % Absolute center frequencies is cross sum of coarse + fine
        % row # index fine freq, col # index coarse freq
        % IMPORTANT: Due to even # sub bands and one centered at 0, the
        % highest and lowest half band share the same band. That bin is
        % assumed to be the highest fine frequency bin.
        freq = repmat(freq,1,length(coarse_freq))+repmat(coarse_freq.',length(freq),1);
        
        % Now merge fine and coarse dimensions together
        spec = permute(spec,[1,3,2]);
        spec = reshape(spec,[],size(spec,3));
        % Do exact same for freq matrix to keep track of frequency bins
        freq = reshape(freq,[],size(freq,3));
        
        % Quantize to fewer bits if not fixed point and QB(2) is set
        if ~DOFLOAT && (QB(2) > 0 && QB(2) < CONFIG{2}.output_nt.WordLength-1)
            spec = quantize(spec/2^(24-QB(2)),numerictype(1,QB(2),0),'Round','Saturate');
        end
        
        % Store output of 2nd stage
        if DOFLOAT
            output{2}.out = spec;
        else
            output{2}.out = spec.removefimath;
        end
        output{2}.fbins = freq;
        
        % Do more bit histograms if flag is set
        if exist('bmax','var')
            output{2}.bithist_inp = fiBitHist(S_IN,bmax);
            output{2}.bithist_filtout = fiBitHist(X,bmax);
            output{2}.bithist_fftout = fiBitHist(spec,bmax);
        end
    end
end

