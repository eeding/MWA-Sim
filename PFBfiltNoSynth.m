function Y = PFBfiltNoSynth(inp_signal,pfbmat,output_nt)
%PFBFILTNOSYNTH Fixed point-Compatible Polyphase filtering with no synthesis.
%   PFBfiltNoSynth(inp_signal,pfbmat,output_nt) - filters inp_signal using
%   the polyphase filter given by pfbmat, and returns output without
%   combining individual phase channels.
%
%   -> inp_signal must be a 1 x N x M matrix, where N is the time dimension
%      and M is the index for parallel signals that can be computed
%      simultaneously. (If single signal, M = 1)
%   -> pfbmat must be a A x B dimensional matrix where A is the number of
%      channels and B is the taps per channel.Coefficients can be complex 
%      or even fi objects. However, they must be compatible with the input 
%      signal type.
%   -> output_nt is the MATLAB numerictype object to which the output
%      values are cast into. Set output_nt = 0 for no casting allowing the
%      output to take on floating values.
%
%   -> Outputs a A x Q x M matrix, where the first dimension indexes the
%      frequency channel, the second dimension indexes the time, and the
%      third dimension indexes the signal (if input has parallel signals).
%      Q = floor(N/A)-B+1;

    threadcount = size(inp_signal,3);
    order = size(pfbmat,1);
    taps = size(pfbmat,2);
    % truncate signal to nearest multiple of filter window size
    r = mod(size(inp_signal,2),order);
    if r ~= 0
        inp_signal = inp_signal(:,1:end-r,:);
    end
    % flip input for convolution
    % important to do this before reshape so the ordering comes out right
    inp_signal = fliplr(inp_signal);
    % reshape (ie. time multiplex) input signal into polyphase form
    inp_signal = reshape(inp_signal,order,[],threadcount);
    
    % pre-generate output matrix
    outputLength = size(inp_signal,2)-taps+1;
    Y = zeros(order,outputLength,threadcount);
    
    % output_nt specifies fi type of output
    % if 0, assume floating point, so input must be cast to float as well
    if output_nt~=0
       Y = fi(Y,output_nt);
       Y.fimath = pfbmat.fimath;
    else
        inp_signal = double(inp_signal);
    end
    % Need a copy of filter per parallel stream
    pfbmat = repmat(pfbmat,1,1,threadcount);
    
    % Do the filtering (shift, multiply and sum)
    % Normally would shift signal into filter from left side
    % To make indexing easier, flip both filter and signal so the
    % signal shifts into the filter from the right side
    inp_signal = fliplr(inp_signal);
    pfbmat = fliplr(pfbmat);
    for n = 1:outputLength
        Y(:,n,:) = sum((pfbmat.*inp_signal(:,n:n+taps-1,:)),2);
    end
end

