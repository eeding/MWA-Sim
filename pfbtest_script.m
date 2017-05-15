clear;load PFBsettings.mat;

% PFBdata{1}.coeff = ones(size(PFBdata{1}.coeff));
% PFBdata{2}.coeff = ones(size(PFBdata{2}.coeff));

% PFBdata{1}.coeff = double(PFBdata{1}.fi_coeff);
% PFBdata{2}.coeff = double(PFBdata{2}.fi_coeff);

%% Global Parameters
% Flags for options
QUANT_INPUT = true;    % Whether or not to quantize input signal
DO_FIXED_PT = true;    % Do simulation with fixed point calculations
DO_FLOAT_PT = true;     % Do simulation with floating point calculations
DO_BIT_HIST = true;    % Produce bit histograms between channelization stages
ADD_NOISE = false;      % Adds IID white Gaussian noise
SNR = 0.0001;           % Ratio of signal power to noise power, only used if ADD_NOISE = true

% 655.36MHz input sampling rate
FS = 128*512*1e4;   
% Time dimension of desired visibility output
D = 1;
% Duration in seconds over which correlations are to be summed
T = 0.001;
% Number of adjacent Fine Channel bins that are to be summed.
% C must be a factor of 128;
C = 4;
% Number of bits for quantization in channelizer
% QB(1) = quantizer after FFB1, QB(2) = quantizer after PFB2;
% Set value to 0 for no quantization
QB = [0,0];
% From FS, D, and T, can calculate length needed for input signal
L = 512*(128*((1e4)*T*D + 11)+7);

% Relative time delay between signal 1 and signal 2 in seconds
% For now assume delay is always ineger multiple of sampling period
DELAY = 1e6/FS;
TAU = ceil(FS*DELAY);  

%% Cosine Test Signals

TIME = 0:1/FS:(L-1)/FS;
coarse_bin = 60;
fine_bin = 37.8;
Fc = (coarse_bin/512+fine_bin/(512*128))*FS;
% cosine amplitude
A = 5;                  

figtitle.cos = ['S(t)=',num2str(A),'Cos(2\pit(',num2str(round(Fc)),'))'];

% Select 24 coarse channels around expected frequency peak
CHSEL = coarse_bin;
CHSEL = max(CHSEL-11,1):1:min(CHSEL+12,256);
CHGAIN = ones(size(CHSEL));

% Build signals
s1 = A*cos(2*pi*Fc*TIME);
s2 = A*cos(2*pi*Fc*(TIME+DELAY));
% Test complex input (debugging)
% s1 = A*(0.5*exp(2*pi*Fc*1i*TIME)+0.5*exp(-2*pi*Fc*1i*TIME));
% s2 = A*exp(2*pi*Fc*1i*(TIME+DELAY));
clear TIME; % Free up some memory

if ADD_NOISE
    n = randn(2,L);
    n = n*sqrt((A*A/2)/SNR);    % Cosine power = A^2/2
    S1 = s1+n(1,:);
    S2 = s2+n(2,:);
else
    S1 = s1;
    S2 = s2;
end
if QUANT_INPUT
    S1 = fi(S1,1,9,0);
    S2 = fi(S2,1,9,0);
    figtitle.cos = sprintf('(Quantized Input): %s',figtitle.cos);
else
    figtitle.cos = sprintf('(Unquantized Input): %s',figtitle.cos);
end

%% Simulate

if DO_FLOAT_PT
    OUTPUT_DATA.floating.sig = {struct(),struct()};
    OUTPUT_DATA.floating.sig{1}.stage = PFBChannelize(FS,S1,PFBdata,CHSEL,CHGAIN,QB,DO_BIT_HIST,1);
    OUTPUT_DATA.floating.sig{2}.stage = PFBChannelize(FS,S2,PFBdata,CHSEL,CHGAIN,QB,DO_BIT_HIST,1);
    OUTPUT_DATA.floating.visibilities = crossMultSum(...
                            OUTPUT_DATA.floating.sig{1}.stage{2}.out,...
                            OUTPUT_DATA.floating.sig{2}.stage{2}.out,...
                            10e3,T,C);
    OUTPUT_DATA.floating.visfbins = mean(reshape(OUTPUT_DATA.floating.sig{1}.stage{2}.fbins,C,[]),1).';
end
if DO_FIXED_PT
    OUTPUT_DATA.fix.sig = {struct(),struct()};
    OUTPUT_DATA.fix.sig{1}.stage = PFBChannelize(FS,S1,PFBdata,CHSEL,CHGAIN,QB,DO_BIT_HIST,0);
    OUTPUT_DATA.fix.sig{2}.stage = PFBChannelize(FS,S2,PFBdata,CHSEL,CHGAIN,QB,DO_BIT_HIST,0);
    OUTPUT_DATA.fix.visibilities = crossMultSum(...
                            OUTPUT_DATA.fix.sig{1}.stage{2}.out,...
                            OUTPUT_DATA.fix.sig{2}.stage{2}.out,...
                            10e3,T,C);
    OUTPUT_DATA.fix.visfbins = mean(reshape(OUTPUT_DATA.fix.sig{1}.stage{2}.fbins,C,[]),1).';                       
end

% save('.\data\dump.mat','OUTPUT_DATA');

%% Compare PFBs

plotlist = {'floating','fix'};
doplot = [DO_FLOAT_PT,DO_FIXED_PT];

inp1_spectrum = abs(fft(double(S1),512));
inp1_spectrum = inp1_spectrum(1:256);
inp1_spectrum = inp1_spectrum(CHSEL);
inp1_spectrum = inp1_spectrum/max(inp1_spectrum);
inp2_spectrum = abs(fft(double(S1),512*128));
inp2_spectrum = inp2_spectrum(1:512*128/2);
inp2_spectrum = inp2_spectrum/max(inp2_spectrum);

w3 = 0:1e4:((length(inp2_spectrum)-1)*1e4);
w4 = 0:100:((length(inp2_spectrum)-1)*1e4-100);
crs_sinc = abs(sinc((w4-Fc)/(FS/512)));
crs_sinc = crs_sinc/crs_sinc(12800*round(Fc/(FS/512))+1);
fine_sinc = abs(sinc((w4-Fc)/(FS/(128*512))));
fine_sinc = fine_sinc/fine_sinc(100*round(Fc/(FS/(512*128)))+1);

for k = 1:2
    if doplot(k)
        w1 = OUTPUT_DATA.(plotlist{k}).sig{1}.stage{1}.fbins;
        w1 = w1(CHSEL);
        pfb1_spectrum = abs(double(OUTPUT_DATA.(plotlist{k}).sig{1}.stage{1}.out));
        pfb1_spectrum = mean(pfb1_spectrum,2);
        pfb1_spectrum = pfb1_spectrum(CHSEL);
        pfb1_spectrum = pfb1_spectrum/max(pfb1_spectrum);

        pfb2_spectrum = abs(double(OUTPUT_DATA.(plotlist{k}).sig{1}.stage{2}.out));
        pfb2_spectrum = mean(pfb2_spectrum,2);
        pfb2_spectrum = pfb2_spectrum/max(pfb2_spectrum);
        w2 = OUTPUT_DATA.(plotlist{k}).sig{1}.stage{2}.fbins;


        fgh = figure(k*100);clf;
        fgh.set('position',[100,100,800,600]);
        subplot(2,1,1);
        hold on;
        stem(w1,inp1_spectrum,'m+','DisplayName','System Input (512 FFT)');
        plot(w1,pfb1_spectrum,'go','DisplayName',sprintf('Coarse PFB (%sPt) Output',plotlist{k}),'linewidth',2);
        plot(w4,crs_sinc,'DisplayName', 'Coarse Sinc');
        hold off;
        legend('show');
        xlabel("Frequency (Hz)");
        ylabel("Magnitude (Rescaled)");

        subplot(2,1,2);
        hold on;
        stem(w3,inp2_spectrum,'r+','DisplayName','System Input (65536 FFT)');
        plot(w2,pfb2_spectrum,'bx','DisplayName',sprintf('Fine PFB (%sPt) Output',plotlist{k}),'linewidth',2);
        plot(w4,fine_sinc,'DisplayName', 'Fine Sinc','Color','black');
        hold off;
        legend('show');
        xlabel("Frequency (Hz)");
        ylabel("Magnitude (Rescaled)");
        suptitle(figtitle.cos);
    end
end

%%
