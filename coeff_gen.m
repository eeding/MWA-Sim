%% Note
%{
This script was used to generate and package the filter coefficients used
in the simulation
There is also a section below that plots the time and frequency response of
each filter.
%}

%%
clear;
% Generate coefficients for Coarse PFB
N = 512;
k = 8;
beta = 5;
bitlen = 12;
h = fir1(N*k-1,1/N,'low',kaiser(N*k,beta))';
h = h./max(h);
coarsePFBcoeff = reshape(h,N,k);

% Quantize to integer
h = (2^(bitlen-1) - 1)*h;

% Fine PFB coefficents are pre-given, load from file.
rawcoeff = csvread('pfb2coeff.csv');
fi_finePFBcoeff = reshape(rawcoeff,128,12);
finePFBcoeff = fi_finePFBcoeff./max(rawcoeff);

% Define the numeric types for PFB coefficients
nt_coarsepfb = numerictype(1,12,0);
nt_finepfb = numerictype(1,24,0);
% Define the fixed-point math logic for each
fm_coarsepfb = fimath('RoundingMethod','Round',...
                          'OverflowAction','Saturate',...
                          'SumMode','SpecifyPrecision',...
                          'SumWordLength',16,...
                          'SumFractionLength',0,...
                          'ProductMode','SpecifyPrecision',...
                          'ProductWordLength',16,...
                          'ProductFractionLength',0);
fm_finepfb = fimath('RoundingMethod','Round',...
                          'OverflowAction','Saturate',...
                          'SumMode','SpecifyPrecision',...
                          'SumWordLength',24,...
                          'SumFractionLength',0,...
                          'ProductMode','SpecifyPrecision',...
                          'ProductWordLength',24,...
                          'ProductFractionLength',0);
% Cast coefficients into fixed point
fi_coarsePFBcoeff = fi(reshape(h,N,k),nt_coarsepfb);
fi_finePFBcoeff = fi(fi_finePFBcoeff,nt_finepfb);
fi_coarsePFBcoeff.fimath = fm_coarsepfb;
fi_finePFBcoeff.fimath = fm_finepfb;

% Define the twiddle coefficients of the fixed point FFTs
fft512twiddle = fi(fi_radix2twiddles(512),1,12);
fft128twiddle = fi(fi_radix2twiddles(128),1,24); 

% Define the output numeric types for the PFB stage (after FFT)
nt_coarsepfb_out = numerictype(1,16,0);
nt_finepfb_out = numerictype(1,24,0);

% Package everything into a list of two structs 
coarse.coeff = coarsePFBcoeff;
coarse.fi_coeff = fi_coarsePFBcoeff;
coarse.twiddle = fft512twiddle;
coarse.output_nt = nt_coarsepfb_out;
fine.coeff = finePFBcoeff;
fine.fi_coeff = fi_finePFBcoeff;
fine.twiddle = fft128twiddle;
fine.output_nt = nt_finepfb_out;
PFBdata = {coarse, fine};
% save('PFBsettings.mat','PFBdata');

%% Plots

figh = figure(1);clf;
figh.set('Position',[100,0,800,800]);
subplot(2,2,1);
plot(1:512*8,coarsePFBcoeff(:),'-o');
axis tight;
xlabel('Index');
ylabel('Amplitude');
title('Floating Precision Filter Impulse Response');
subplot(2,2,2);
spc = abs(fft(coarsePFBcoeff(:)));
spc = spc(1:end/64)/max(spc);
w = linspace(0,0.5/32,length(spc));
plot(w,10*log10(spc),'-o');
axis tight;
xlabel('Normalized Frequency (Units of F_{samp})');
ylabel('Normalized Amplitude (dB)');
title('Floating Precision Filter Frequency Response');
subplot(2,2,3);
plot(1:512*8,fi_coarsePFBcoeff(:),'-o');
axis tight;
xlabel('Index');
ylabel('Amplitude');
title('Integer Filter Impulse Response');
subplot(2,2,4);
spc = abs(fft(double(fi_coarsePFBcoeff(:))));
spc = spc(1:end/64)/max(spc);
w = linspace(0,0.5/32,length(spc));
plot(w,10*log10(spc),'-o');
axis tight;
xlabel('Normalized Frequency (Units of F_{samp})');
ylabel('Normalized Amplitude (dB)');
title('Integer Filter Frequency Response');

suptitle('Coarse PFB Prototype Filter');



figh = figure(2);clf;
figh.set('Position',[100,0,800,800]);
subplot(2,2,1);
plot(1:128*12,finePFBcoeff(:),'-o');
axis tight;
xlabel('Index');
ylabel('Amplitude');
title('Floating Precision Filter Impulse Response');
subplot(2,2,2);
spc = abs(fft(finePFBcoeff(:)));
spc = spc(1:end/64)/max(spc);
w = linspace(0,0.5/32,length(spc));
plot(w,10*log10(spc),'-o');
axis tight;
xlabel('Normalized Frequency (Units of F_{samp})');
ylabel('Normalized Amplitude (dB)');
title('Floating Precision Filter Frequency Response');
subplot(2,2,3);
plot(1:128*12,fi_finePFBcoeff(:),'-o');
axis tight;
xlabel('Index');
ylabel('Amplitude');
title('Integer Filter Impulse Response');
subplot(2,2,4);
spc = abs(fft(double(fi_finePFBcoeff(:))));
spc = spc(1:end/64)/max(spc);
w = linspace(0,0.5/32,length(spc));
plot(w,10*log10(spc),'-o');
axis tight;
xlabel('Normalized Frequency (Units of F_{samp})');
ylabel('Normalized Amplitude (dB)');
title('Integer Filter Frequency Response');

suptitle('Fine PFB Prototype Filter');



figh = figure(3);clf;
sq1 = ones(1,512*8);
sq2 = ones(1,128*12);
figh.set('Position',[100,0,800,800]);
subplot(2,2,1);
plot(1:512*8,sq1,'-o');
axis tight;
xlabel('Index');
ylabel('Amplitude');
title('Square Filter 1 Impulse Response');
subplot(2,2,2);
spc = abs(fft(sq1));
spc = spc(1:end)/max(spc);
w = linspace(0,0.5/32,length(spc));
plot(w,spc,'-o');
axis tight;
xlabel('Normalized Frequency (Units of F_{samp})');
ylabel('Normalized Amplitude');
title('Square Filter 1 Frequency Response');
subplot(2,2,3);
plot(1:128*12,sq2,'-o');
axis tight;
xlabel('Index');
ylabel('Amplitude');
title('Square Filter 2 Impulse Response');
subplot(2,2,4);
spc = abs(fft(double(sq2)));
spc = spc(1:end)/max(spc);
w = linspace(0,0.5/32,length(spc));
plot(w,spc,'-o');
axis tight;
xlabel('Normalized Frequency (Units of F_{samp})');
ylabel('Normalized Amplitude');
title('Square Filter 2 Frequency Response');

suptitle('Square Prototype Filter');
