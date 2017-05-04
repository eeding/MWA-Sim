function [output] = crossMultSum(sig1,sig2,FS,T,C)
%CROSSMULTSUM simulates the accumulation and summation stage of the MWA DSP
%             chain.
%
%   crossMultSum(sig1,sig2,FS,T,C) - pointwise multiplies sig1 and sig2
%   which are both sampled at FS Hz. The product is sum over T seconds and
%   over every C adjacent channels.
%
%   -> sig1 - M x N matrix where every row corresponds to a partiular
%             signal channel and every column is a time index.
%   -> sig2 - Another matrix with the same specification as sig1. The
%             pointwise product of sig2 and sig1 will be computed.
%   -> FS - The sampling frequency in Hz of all the input channels in both
%           signal sets.
%   -> T - The duration in seconds over which the product of sig1 and sig2
%          should be summed.
%   -> C - The number of adjacent channels over which the product of sig1
%          and sig2 should be summed. C must divide 128 since the
%          simulation assumes M = 3072 = 24 sets of 128 fine channels.
    
    if ~isequal(size(sig1),size(sig2))
        error('ERROR: Signal dimensions unequal!');
    else
        if (~isequal(size(size(sig1)),[1,2]))
            error('ERROR: Signals must be 2D (row = channel, column = time');
        end
    end
    
    % Make sure signals are floating point before computing
    % Conjugate one of them in case complex value
    sig1 = double(sig1);
    sig2 = double(sig2);
    
    % Perform pointwise multiplication
    xprod = sig1.*sig2;
    % Rotate channel dimension "out of the way" to dimension 3
    % in peparation for parallel summation
    xprod = permute(xprod,[3,2,1]);
    
    % Calculate how many time samples to sum over
    sumsize = T*FS;
    if 0 ~= mod(sumsize,1)
        disp('WARNING: Selected sum period is not multiple of sampling period.');
        sumsize = ceil(T*FS);
        disp(strcat('Using T = ',num2str(sumsize/FS),'s instead.'));
    end
    r = mod(size(xprod,2),sumsize);
    if (r ~= 0)
        disp('WARNING: Signal duration does not evenly divide sum period.');
        disp(strcat('Last data point in output will be over T = ',num2str(r/FS),'s'));
        % pad with 0 which will not affect sum
        xprod = [xprod,zeros(1,sumsize-r,size(xprod,3))];
    end
    % reshape so elements to be summed are in same column
    output = reshape(xprod,sumsize,[],size(xprod,3));
    output = sum(output,1);

    % Optionally sum over adjacent C channels
    if (C > 0) && (mod(128,C) == 0)
        % Rotate back channel dimension, and rotate time dimension away
        output = permute(output,[3,1,2]);
        % Reshape to sum over consecutive C channels
        output = reshape(output,C,[],size(output,3));
        output = sum(output,1);
        % Final output has the channels as rows and the time as columns
        output = permute(output,[2,3,1]);
    else
        % otherwise rotate just back channel dimension
        output = permute(output,[3,2,1]);
    end
        
end

