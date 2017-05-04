function [S,f] = fi_radix2fft(X,W)
%FI_RADIX2FFT  Radix-2 Fixed point FFT
%   fi_radix2fft(X, W) - computes the fixed point radix-2 FFT across the
%   first dimension (rows) of the input matrix X using the twiddle factors
%   given by W.
%
%   ->  X and W is assumed to contain fi objects.
%   ->  The number of rows (N) in X must be a power of 2.
%   ->  W must be the twiddle factors corresponding to an FFT of order N. 
%   ->  Use fi_radix2twiddles.m to generate floating point twiddles.
%       W = fi_radix2twiddles(N)
%       The floating point W can then be converted into fi objects.
%
%   -> Outputs a tuple S,f.
%      S is a matrix with the same dimension as X, but now containing the
%      FFT calculated along the first dimension.
%      f gives the normalized frequency value of each frequency bin for
%      the FFT. The jth element of f is the bin for the jth index along
%      dim 1 of S.
% 
%   This function was based on the implementation given here:
%   https://www.mathworks.com/help/fixedpoint/ug/convert-fast-fourier-transform-fft-to-fixed-point.html
   

    S = X;
    for p = 1:size(X,3)
        S(:,:,p) = bitrevorder(X(:,:,p));
    end
   
    % pre-generate index variables
    n = size(S,1);
    t = log2(n);
    LL = cast(2.^(1:t),'like',int32([]));
    rr = cast(n./LL,'like',int32([]));
    LL2 = cast(LL./2,'like',int32([]));
    % Compute fft in parallel over all columns and threads
    for q=1:t
        L = LL(q); r = rr(q); L2 = LL2(q);
        for k=0:L:(L*(r-1))
            for j=1:L2
                % Compute butterflies
                temp = W(L2-1+j)*S(k+j+L2,:,:);
                % bitshift normalizes power by input length
                S(k+j+L2,:,:) = bitsra(S(k+j,:,:)-temp,1);
                S(k+j,:,:)= bitsra(S(k+j,:,:)+temp,1);
            end
        end
    end
    % rearrange so frequency to go from 0 to 1 (normalized)
    S = [S(1,:,:);flipud(S(n/2+1:end,:,:)); flipud(S(2:n/2,:,:))];
    % label frequencies above 0.5 with negative equivalent
    f = [(0:1/n:0.5),((-0.5+1/n):1/n:-1/n)]';
end