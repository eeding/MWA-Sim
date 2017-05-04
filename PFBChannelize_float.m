function [output] = PFBChannelize_float(FS,S_IN,CONFIG,CHSEL,CHGAIN)
%PFBCHANNELIZE_FLOAT uses PFBFilterNoSynth to perform the two stage PFB 
%                    filtering in the floating point MWA DSP simulation.
%
%   This function operates in the same way as PFBCHANNELIZE, except it
%   operates on floating point input and performs no quantization and does
%   not have any quantization logic in the filtering process. It makes use
%   of MATLAB's built in FFT function instead of fi_radix2fft.

    output = {struct(), struct()};
    n1 = size(CONFIG{1}.coeff,1); % this should be 512
    n2 = size(CONFIG{2}.coeff,1); % this should be 128
    if (~isequal(size(size(CHSEL)),[1,2]))
        error('ERROR: CHSEL must be row vector or 0');
    else
        CHSEL = mod(abs(CHSEL-1),ceil(n1/2))+1;
    end
    if ~isequal(size(CHGAIN),size(CHSEL))
        disp('CHGAIN dimension does not match CHSEL dimension. Assuming unity gain.');
        CHGAIN = ones(size(CHSEL));
    end
    
    X = PFBfiltNoSynth(S_IN,CONFIG{1}.coeff,0);
    if size(X,2) < 1
       error('Insufficient data.');
    end
    spec = fft(X,n1,1)./n1;
    freq = -0.5:1/n1:0;
    freq = [0:1/n1:0.5,freq(2:end-1)]';
    % store stage 1 data
    output{1}.out = spec(1:ceil(n1/2),:,:);
    output{1}.fbins = FS.*freq(1:ceil(n1/2),:,:);
    
    % Check if there's enough datapoints to continue
    if(size(output{1}.out,2) < length(CONFIG{2}.coeff(:)))
        disp("Not enough samples to proceed to stage 2!");
        return;
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
        X = PFBfiltNoSynth(X,CONFIG{2}.coeff,0);
        spec = fft(X,n2,1)./n2;
        freq = -0.5:1/n2:0;
        freq = [0:1/n2:0.5,freq(2:end-1)]';
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
        % Store output of 2nd stage
        output{2}.out = spec;
        output{2}.fbins = freq;
    end
end

