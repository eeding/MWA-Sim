clear;load sim_settings.mat;

%% Global Parameters
% Flags for options
QUANT_INPUT = false;    % Whether or not to quantize input signal
DO_FIXED_PT = false;    % Do simulation with fixed point calculations
DO_FLOAT_PT = true;     % Do simulation with floating point calculations
DO_BIT_HIST = false;    % Produce bit histograms between channelization stages
ADD_NOISE = false;      % Adds IID white Gaussian noise
SNR = 0.0001;           % Ratio of signal power to noise power, only used if ADD_NOISE = true
PLOT_DIM = 3;           % How many dimensional plot? (2D or 3D are choices)

% 655.36MHz input sampling rate
FS = 128*512*1e4;   
% Time dimension of desired visibility output
D = 10;
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
Fc = FS/12;         % cosine frequency
A = 5;              % cosine amplitude

% Select 24 coarse channels around expected frequency peak
CHSEL = floor(Fc/(FS/512));
CHSEL = max(CHSEL-11,1):1:min(CHSEL+12,256);
CHGAIN = ones(size(CHSEL));

% Build signals
s1 = A*cos(2*pi*Fc*TIME);
s2 = A*cos(2*pi*Fc*(TIME+DELAY));
clear TIME; % Free up some memory

% Pre-generate title for plotting later
fig_title = ['S(t)=',num2str(A),'Cos(2\pit(',num2str(round(Fc)),'))'];

if ADD_NOISE
    n = randn(2,L);
    n = n*sqrt((A*A/2)/SNR);    % Cosine power = A^2/2
    S1 = s1+n(1,:);
    S2 = s2+n(2,:);
    fig_title =[fig_title,' | With Noise (SNR=',num2str(SNR),')'];
else
    S1 = s1;
    S2 = s2;
end
if QUANT_INPUT
    S1 = fi(S1,1,9,0);
    S2 = fi(S2,1,9,0);
    fig_title = ['[Quantized] ', fig_title];
end

fig_title = strcat(fig_title,' | \tau = ',num2str(DELAY*1e3),'ms |');
%% Random Guassian Input Signals

% Input signal will be bandpassed to this interval (Hz)
% Set to 0 for no input filtering
BAND_LIM = [80e6,300e6];    % 80MHz to 300MHz
BAND_LIM = BAND_LIM/(FS/2);
% Overall signal gain
A = 4;

% Pick 24 coarse channel around a random bin
center = randi([50,200]);
CHSEL = max(center-11,1):1:min(center+12,256);
CHGAIN = ones(size(CHSEL));

% Build Signals
sigma = sqrt(SNR);
s1 = randn(1,L)*sigma;
s2 = [zeros(1,TAU),s1(1:end-TAU)];

% Pre-generate title for plotting later
fig_title = strcat('S(t)~N(\mu=0,\sigma=',num2str(A*sigma),')');

if ADD_NOISE
    n = randn(2,L);
    S1 = s1 + n(1,:);
    S2 = s2 + n(2,:);
    fig_title =[fig_title,' | With Noise (SNR=',num2str(SNR),')'];
else
    S1 = s1;
    S2 = s2;
end
% Overall gain applied after adding noise
S1 = A*S1;
S2 = A*S2;

% Band limit input signals
if boolean(BAND_LIM)
    filtspecs = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
    BAND_LIM(1),BAND_LIM(1)*1.05,BAND_LIM(2)*0.95,BAND_LIM(2),60,1,60);
    bpf = design(filtspecs,'equiripple');
    S1 = filter(bpf,S1);
    S2 = filter(bpf,S2);
end

if QUANT_INPUT
    S1 = fi(S1,1,9,0);
    S2 = fi(S2,1,9,0);
    fig_title = ['[Quantized] ', fig_title];
end

fig_title = strcat(fig_title,' | \tau = ',num2str(DELAY*1e3),'ms |');
%% Simulate

if DO_FLOAT_PT
    OUTPUT_DATA.flt.sig = {struct(),struct()};
    OUTPUT_DATA.flt.sig{1}.stage = PFBChannelize_float(FS,double(S1),PFBdata,CHSEL,CHGAIN);
    OUTPUT_DATA.flt.sig{2}.stage = PFBChannelize_float(FS,double(S2),PFBdata,CHSEL,CHGAIN);
    OUTPUT_DATA.flt.visibilities = crossMultSum(...
                            OUTPUT_DATA.flt.sig{1}.stage{2}.out,...
                            OUTPUT_DATA.flt.sig{2}.stage{2}.out,...
                            10e3,T,C);
    OUTPUT_DATA.flt.visfbins = mean(reshape(OUTPUT_DATA.flt.sig{1}.stage{2}.fbins,C,[]),1)';
end
if DO_FIXED_PT
    OUTPUT_DATA.fix.sig = {struct(),struct()};
    OUTPUT_DATA.fix.sig{1}.stage = PFBChannelize(FS,S1,PFBdata,CHSEL,CHGAIN,QB,DO_BIT_HIST);
    OUTPUT_DATA.fix.sig{2}.stage = PFBChannelize(FS,S2,PFBdata,CHSEL,CHGAIN,QB,DO_BIT_HIST);
    OUTPUT_DATA.fix.visibilities = crossMultSum(...
                            OUTPUT_DATA.fix.sig{1}.stage{2}.out,...
                            OUTPUT_DATA.fix.sig{2}.stage{2}.out,...
                            10e3,T,C);
    OUTPUT_DATA.fix.visfbins = mean(reshape(OUTPUT_DATA.fix.sig{1}.stage{2}.fbins,C,[]),1)';                       
end

% save('dump.mat','OUTPUT_DATA');

%% Plots

fig_titles{1} = strcat(fig_title,' Floating Precision');
fig_titles{2} = strcat(fig_title,' Fixed Precision');
fig_titles{3} = ['Visibilities (T=',num2str(T),'s, C=',num2str(C),'): ',fig_titles{1}];
fig_titles{4} = ['Visibilities (T=',num2str(T),'s, C=',num2str(C),'): ',fig_titles{2}];

selection = {'flt','fix'};
prec = [DO_FLOAT_PT,DO_FIXED_PT];
for k = 1:2
    if(prec(k))
        currdata = OUTPUT_DATA.(selection{k});
        fighandle = figure(k);clf;
        fighandle.set('Position',[100,50,1000,700]);
        for stg = 1:2
            Xscale = 512/FS; 
            for n = 1:2 
                subplot(2,2,(stg-1)*2+n);
                % plotting dB
                % Adding 1 offset to zero out everything below 0 dB
                Z = 10*log10(abs(double(currdata.sig{n}.stage{stg}.out)));
                X = currdata.sig{n}.stage{stg}.fbins;
                Y = (1:size(Z,2))*Xscale;
                surf(X',Y',Z','EdgeColor','None','Facecolor','Interp');
                view(PLOT_DIM);
                colormap('Jet');colorbar;axis tight;
                xlabel(colorbar(),'Amplitude (dB)');
                title(['PFB ',num2str(stg),' output for Signal ',num2str(n)]);
                xlabel('Frequency (Hz)');
                ylabel('Time (s)');
            end
            Xscale = Xscale*128;
        end
        suptitle(fig_titles{k});

        % Plot Visibilies
        fighandle = figure(10*k);clf;
        fighandle.set('Position',[100,50,1000,600]);
        Z = 10*log10(abs(currdata.visibilities));
        if size(Z,2) < 2
            % if only 1 visibility time sample, duplicate so can 3D plot
            Z = repmat(Z,1,2);
        end
        X = currdata.visfbins;
        Y = (1:size(Z,2))*T;
        surf(X',Y',Z','EdgeColor','None','Facecolor','Interp');
        view(PLOT_DIM);
        colormap('Jet');colorbar;axis tight;
        xlabel(colorbar(),'Amplitude (dB)');
        suptitle(fig_titles{k+2});
        xlabel('Frequency Bin Centers (Hz)');
        ylabel('Time (s)');
    end
end
