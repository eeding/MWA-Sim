%% Set Parameters

clear;
% Set true if debugging with cosine. check cosine code below tweak what
% cosine is being generated as the input
% If set false, will use Gaussian input with given SNR and noise sigma = 4
SIM_PARAM.USE_COS_INPUT = true;

% set true if want to save output data struct to file
SAVE_OUTPUT = false;

% Flags for options
SIM_PARAM.QUANT_INPUT = false;    % Whether or not to quantize input signal
SIM_PARAM.DO_FIXED_PT = false;    % Do simulation with fixed point calculations
SIM_PARAM.DO_FLOAT_PT = true;     % Do simulation with floating point calculations
SIM_PARAM.DO_BIT_HIST = false;    % Produce bit histograms between channelization stages
SIM_PARAM.ADD_NOISE = false;      % Adds IID white Gaussian noise
SIM_PARAM.INPUT_SNR = 0.0001;     % Ratio of signal power to noise power, only used if ADD_NOISE = true

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

% Coarse channel selected after 1st stage PFB. Must be in range [1,256]
% 1 is baseband, 256 is highest band.
% Default is pick 24 coarse channel around a random bin
center = randi([50,200]);
CHSEL = max(center-11,1):1:min(center+12,256);
% Gain for each coarse channel, default is 1
CHGAIN = ones(size(CHSEL));

% Relative time delay between signal 1 and signal 2 in seconds
% For now assume delay is always ineger multiple of sampling period
DELAY = 1e6/FS;
TAU = ceil(FS*DELAY);  

if SAVE_OUTPUT
    SIM_PARAMS.VIS_SUM_PERIOD = T;
    SIM_PARAMS.VIS_CHAN_BW = C*1e4;
    SIM_PARAMS.PFB1_QUANT_BITSIZE  = QB(1);
    SIM_PARAMS.PFB2_QUANT_BITSIZE  = QB(2);
    SIM_PARAMS.INP_REL_DELAY = DELAY;
    SAVEFILE_PATH = strcat('MWA_SIM_',datestr(now,'mm-dd-yyyy(HH-MM)'));
    SAVEFILE_PATH = strcat('./data/',SAVEFILE_PATH,'.mat');
end


%% Run Simulation

%%%======================================================%%%
%%%=================== INPUT CREATION ===================%%%
%%%======================================================%%%
if SIM_PARAM.USE_COS_INPUT
    % Generate cosine test signal
    
    TIME = 0:1/FS:(L-1)/FS;
    Fc = FS/12;         % cosine frequency
    A = 5;              % cosine amplitude

    % Select 24 coarse channels around expected frequency peak
    CHSEL = floor(Fc/(FS/512));
    CHSEL = max(CHSEL-11,1):1:min(CHSEL+12,256);
    CHGAIN = ones(size(CHSEL));

    % Build signals
    SIM_INPUT.OF_INTEREST{1} = A*cos(2*pi*Fc*TIME);
    SIM_INPUT.OF_INTEREST{2} = A*cos(2*pi*Fc*(TIME+DELAY));
    clear TIME; % Free up some memory

    % Pre-generate title for plotting later
    fig_title = ['S(t)=',num2str(A),'Cos(2\pit(',num2str(round(Fc)),'))'];

    if SIM_PARAM.ADD_NOISE
        n = randn(2,L);
        n = n*sqrt((A*A/2)/SIM_PARAM.INPUT_SNR);    % Cosine power = A^2/2
        SIM_INPUT.RAW{1} = SIM_INPUT.OF_INTEREST{1}+n(1,:);
        SIM_INPUT.RAW{2} = SIM_INPUT.OF_INTEREST{2}+n(2,:);
        fig_title =[fig_title,' | With Noise (SNR=',num2str(SIM_PARAM.INPUT_SNR),')'];
    else
        SIM_INPUT.RAW{1} = SIM_INPUT.OF_INTEREST{1};
        SIM_INPUT.RAW{2} = SIM_INPUT.OF_INTEREST{2};
    end
    if SIM_PARAM.QUANT_INPUT
        SIM_INPUT.RAW{1} = fi(SIM_INPUT.RAW{1},1,9,0);
        SIM_INPUT.RAW{2} = fi(SIM_INPUT.RAW{2},1,9,0);
        fig_title = ['[Quantized] ', fig_title];
    end

    fig_title = strcat(fig_title,' | \tau = ',num2str(DELAY*1e3),'ms |');
else
    % Generate Random Guassian input signal
    % Input signal will be bandpassed to this interval (Hz)
    % Set to 0 for no input filtering
    BAND_LIM = [80e6,300e6];    % 80MHz to 300MHz
    BAND_LIM = BAND_LIM/(FS/2);
    % Overall signal gain
    A = 4;

    % Build Signals
    sigma = sqrt(SIM_PARAM.INPUT_SNR);
    SIM_INPUT.OF_INTEREST{1} = randn(1,L)*sigma;
    SIM_INPUT.OF_INTEREST{2} = [zeros(1,TAU),SIM_INPUT.OF_INTEREST{1}(1:end-TAU)];

    % Pre-generate title for plotting later
    fig_title = strcat('S(t)~N(\mu=0,\sigma=',num2str(A*sigma),')');

    if SIM_PARAM.ADD_NOISE
        n = randn(2,L);
        SIM_INPUT.RAW{1} = SIM_INPUT.OF_INTEREST{1} + n(1,:);
        SIM_INPUT.RAW{2} = SIM_INPUT.OF_INTEREST{2} + n(2,:);
        fig_title =[fig_title,' | With Noise (SNR=',num2str(SIM_PARAM.INPUT_SNR),')'];
    else
        SIM_INPUT.RAW{1} = SIM_INPUT.OF_INTEREST{1};
        SIM_INPUT.RAW{2} = SIM_INPUT.OF_INTEREST{2};
    end
    % Overall gain applied after adding noise
    SIM_INPUT.RAW{1} = A*SIM_INPUT.RAW{1};
    SIM_INPUT.RAW{2} = A*SIM_INPUT.RAW{2};

    % Band limit input signals
    if boolean(BAND_LIM)
        filtspecs = fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
        BAND_LIM(1),BAND_LIM(1)*1.05,BAND_LIM(2)*0.95,BAND_LIM(2),60,1,60);
        bpf = design(filtspecs,'equiripple');
        SIM_INPUT.RAW{1} = filter(bpf,SIM_INPUT.RAW{1});
        SIM_INPUT.RAW{2} = filter(bpf,SIM_INPUT.RAW{2});
    end

    if SIM_PARAM.QUANT_INPUT
        SIM_INPUT.RAW{1} = fi(SIM_INPUT.RAW{1},1,9,0);
        SIM_INPUT.RAW{2} = fi(SIM_INPUT.RAW{2},1,9,0);
        fig_title = ['[Quantized] ', fig_title];
    end

    fig_title = strcat(fig_title,' | \tau = ',num2str(DELAY*1e3),'ms |');
end

fig_titles{1} = strcat(fig_title,' Floating Precision');
fig_titles{2} = strcat(fig_title,' Fixed Precision');
fig_titles{3} = ['Visibilities (T=',num2str(T),'s, C=',num2str(C),'): ',fig_titles{1}];
fig_titles{4} = ['Visibilities (T=',num2str(T),'s, C=',num2str(C),'): ',fig_titles{2}];

%%%==================================================%%%
%%%=================== SIMULATION ===================%%%
%%%==================================================%%%
    
if ~exist('PFBsettings','var')
    load PFBsettings.mat;
end
if SIM_PARAM.DO_FLOAT_PT
    SIM_OUTPUT.floatsim.sig = {struct(),struct()};
    SIM_OUTPUT.floatsim.sig{1}.stage = PFBChannelize_float(FS,double(SIM_INPUT.RAW{1}),PFBdata,CHSEL,CHGAIN);
    SIM_OUTPUT.floatsim.sig{2}.stage = PFBChannelize_float(FS,double(SIM_INPUT.RAW{2}),PFBdata,CHSEL,CHGAIN);
    SIM_OUTPUT.floatsim.visibilities = crossMultSum(...
                            SIM_OUTPUT.floatsim.sig{1}.stage{2}.out,...
                            SIM_OUTPUT.floatsim.sig{2}.stage{2}.out,...
                            10e3,T,C);
    SIM_OUTPUT.floatsim.visfbins = mean(reshape(SIM_OUTPUT.floatsim.sig{1}.stage{2}.fbins,C,[]),1)';
end
if SIM_PARAM.DO_FIXED_PT
    SIM_OUTPUT.intsim.sig = {struct(),struct()};
    SIM_OUTPUT.intsim.sig{1}.stage = PFBChannelize(FS,SIM_INPUT.RAW{1},PFBdata,CHSEL,CHGAIN,QB,SIM_PARAM.DO_BIT_HIST);
    SIM_OUTPUT.intsim.sig{2}.stage = PFBChannelize(FS,SIM_INPUT.RAW{2},PFBdata,CHSEL,CHGAIN,QB,SIM_PARAM.DO_BIT_HIST);
    SIM_OUTPUT.intsim.visibilities = crossMultSum(...
                            SIM_OUTPUT.intsim.sig{1}.stage{2}.out,...
                            SIM_OUTPUT.intsim.sig{2}.stage{2}.out,...
                            10e3,T,C);
    SIM_OUTPUT.intsim.visfbins = mean(reshape(SIM_OUTPUT.intsim.sig{1}.stage{2}.fbins,C,[]),1)';                       
end

disp("=========SIMULATION COMPLETE=========")
if SAVE_OUTPUT
    disp("Saving Data...")
    SIM_PARAMS.COARSE_CHAN_SEL = CHSEL;
    SIM_PARAMS.COARSE_CHAN_GAIN = CHGAIN;
    save(SAVEFILE_PATH,'SIM_OUTPUT','SIM_PARAMS','SIM_INPUT','fig_titles');
    disp("Data saved!")
end

%% Make Plots

% How many dimensional plot? (2D or 3D are choices)
PLOT_DIM = 3;     

selection = {'floatsim','intsim'};
prec = [SIM_PARAM.DO_FLOAT_PT,SIM_PARAM.DO_FIXED_PT];
for k = 1:2
    if(prec(k))
        currdata = SIM_OUTPUT.(selection{k});
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
