function [h] = fiBitHist(X,B)
%FIBITHIST returns a 2 x B vector containing the normalized histogram of
%   the binary digits in matrix X. Row 1 is the histogram for 1s, row 2 is
%   the histogram for 0s.
%   fibithist(X,B)- counts the occurence of every binary digit in the data
%   provided in matrix X, up to the Bth MSB. The data in X is assumed to be
%   fi objects which store binary integers as two's complement. The raw
%   bits will be converted to binary interpretation before being
%   histogrammed.
%   -> X must be a matrix containing fix point (fi) objects, each with the
%      same bit length.
%   -> B must be a single postive integer indicating the maximum bits to
%      count up to.

    X = X(:);
    X = [real(X);imag(X)];
    L = length(X);
    h = zeros(2,B);
    for b = 1:B
        try
            h(1,b) = sum(double(bitget(X,b)))/L;
            h(2,b) = 1-h(1,b);
        catch
            break;
        end
    end
end

