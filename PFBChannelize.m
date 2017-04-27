function [output] = PFBChannelize(FS,S_IN,CONFIG,CHSEL,CHGAIN,QB,bithist)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    output = {struct(),struct()};

    n1 = size(CONFIG{1}.fi_coeff,1); % this should be 512
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
   
    % First stage straightforward filtering
    X = PFBfiltNoSynth(S_IN,CONFIG{1}.fi_coeff,CONFIG{1}.output_nt);
    if size(X,2) < 1
        error('Insufficient data.');
    end
    [spec,freq] = fi_radix2fft(X,CONFIG{1}.twiddle);

    % Quantize to fewer bits if QB(1) is set
    if (QB(1) > 0 && QB(1) < 15)
        S_IN = quantize(S_IN/2^(16-QB(1)),numerictype(1,QB(1),0),'Round','Saturate');
    end
    
    % Store results for first stage
    output{1}.out = spec.removefimath;
    output{1}.out = output{1}.out(1:ceil(n1/2),:,:);
    output{1}.fbins = FS.*freq(1:ceil(n1/2),:,:);
    
    % Calculate the bit histograms for stage 1
    % Going up to 24 bits
    bmax = 24;
    if bithist
        output{1}.bithist_inp = fiBitHist(S_IN,bmax);
        output{1}.bithist_filtout = fiBitHist(X,bmax);
        output{1}.bithist_spec = fiBitHist(spec,bmax);
    end
    
    % Check if there's enough datapoints to continue to stage 2
    if(size(output{1}.out,2) < length(CONFIG{2}.coeff(:)))
        disp("Not enough samples to proceed to stage 2!");
        return
    else
        % Select only a few coarse channels for next stage
        S_IN = output{1}.out(CHSEL,:,:);
        coarse_freq = output{1}.fbins(CHSEL);
        % Scale each channel by appropriate gain
        S_IN = S_IN.*repmat(CHGAIN',1,size(S_IN,2),size(S_IN,3));
        
        % Move coarse channel index to dimension 3 so they will be acted on
        % in parallel for second stage
        X = permute(S_IN,[3,2,1]);
        % Perform second stage filtering
        X = PFBfiltNoSynth(X,CONFIG{2}.fi_coeff,CONFIG{2}.output_nt);
        [spec,freq] = fi_radix2fft(X,CONFIG{2}.twiddle);
        % 2nd stage frequency bin is relative to each coarse bin
        freq = (FS/n1).*freq;
        % Reorder to make sure frequencies come out in order
        spec = [spec(freq<0,:,:);spec(freq>=0,:,:)];
        freq = [freq(freq<0);freq(freq>=0)];
        
        % Absolute center frequencies is cross sum of coarse + fine
        % row # index fine freq, col # index coarse freq
        freq = repmat(freq,1,length(coarse_freq))+repmat(coarse_freq',length(freq),1);
        % Now merge fine and coarse dimensions together
        spec = permute(spec,[1,3,2]);
        spec = reshape(spec,[],size(spec,3));
        % Do exact same for freq matrix to keep track of frequency bins
        freq = reshape(freq,[],size(freq,3));
        
        % Quantize to fewer bits if QB(2) is set
        if (QB(2) > 0 && QB(2) < 24)
            spec = quantize(spec/2^(24-QB(2)),numerictype(1,QB(2),0),'Round','Saturate');
        end
        
        % Store output of 2nd stage
        output{2}.out = spec;
        output{2}.fbins = freq;
        
        % Do more bit histograms if flag is set
        if bithist
            output{2}.inp = fiBitHist(S_IN,bmax);
            output{2}.bithist_filtout = fiBitHist(X,bmax);
            output{2}.bithist_spec = fiBitHist(spec,bmax);
        end
    end
end

