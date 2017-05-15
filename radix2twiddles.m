function W = fi_radix2twiddles(n) 
%FI_RADIX2TWIDDLES  Twiddle factors for radix-2 FFT.
%   FI_RADIX2TWIDDLES(N) - computes the length N-1 vector W of
%   twiddle factors to be used in the FI_M_RADIX2FFT function.
%
%   This function was based on the implementation given here:
%   https://www.mathworks.com/help/fixedpoint/ug/convert-fast-fourier-transform-fft-to-fixed-point.html


t = log2(n);
if floor(t) ~= t
  error('N must be an exact power of two.');
end

W = zeros(n-1,1);
k=1;
L=2;
% Equation 1.4.11, p. 34
while L<=n
  theta = 2*pi/L;
  % Algorithm 1.4.1, p. 23
  for j=0:(L/2 - 1)
    W(k) = complex( cos(j*theta), -sin(j*theta) );
    k = k + 1;
  end
  L = L*2;
end