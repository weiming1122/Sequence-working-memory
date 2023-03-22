clc;
clear;
close all;

maindir = pwd;                    % keep main path
cd E:\lwm\fieldtrip-master        % set up the path of fieldtrip
addpath(pwd)

ft_defaults

cd(maindir)                       % return to main

%% setup directories
root = pwd; out = 'data_analysis';
dRoot = [root(1:end-length(out)),'raw_data/'];

%% inspect raw data in different sets of electrodes
cfg = [];
cfg.dataset = [dRoot,'direction1.edf'];
cfg.continuous = 'yes';
cfg.viewmode   = 'butterfly'; % try 'butterfly','vertical'

%% inspect raw data in *A* electrodes (Central)
cfg.channel = {'*A*','-EDF*'};
ft_databrowser(cfg);                      % electrodes:

%% inspect raw data in *B* electrodes (Frontal)
cfg.channel = {'*B*','-*BP*'};
ft_databrowser(cfg);                     % electrodes: POL B30

%% inspect raw data in *C* electrodes (Temporal)
cfg.channel = {'*C*','-*DC*'};
ft_databrowser(cfg);                     % electrodes: 

%% inspect raw data in *D* electrodes (Wernicke)
cfg.channel = {'*D*','-*DC*','-EDF*'};
ft_databrowser(cfg);                     % electrodes:

%% inspect raw data in *E* electrodes (Broca)
cfg.channel = {'POL E*','-POL E','-POL EKG*','-POL EMG*'};
ft_databrowser(cfg);                    % electrodes: POL E23

%% inspect raw data in *DC* electrodes (Trigger)
cfg.channel = {'POL DC*'};
ft_databrowser(cfg);                    % electrodes: 

%% inspect raw data in *EKG* electrodes (Heart)
cfg.channel = {'POL EKG1','POL EKG2'};
ft_databrowser(cfg);                    % electrodes: POL EKG1

%% inspect raw data in *EMG* electrodes (Muscle)
cfg.channel = {'POL EMG1','POL EMG2'};
ft_databrowser(cfg);                    % electrodes: 

%% inspect raw data in *BP* electrodes (?)
cfg.channel = {'POL BP*'};
ft_databrowser(cfg);                    % electrodes: 

%% inspect raw data in *EDF* electrodes (?)
cfg.channel = {'EDF*'};
ft_databrowser(cfg);                    % electrodes: 

%% import trigger data based on trigger channel
cfg = [];
cfg.dataset = [dRoot,'direction1.edf'];
cfg.continuous = 'yes';
cfg.channel = {'*DC10','*DC11','*DC12','*DC13','*DC14'};
dataTrigger = ft_preprocessing(cfg);

% visually inspect the data
cfg = [];
cfg.viewmode = 'vertical';
cfg.ylim = 'maxmin';
ft_databrowser(cfg, dataTrigger);

% extract stimuli01 marker
tm = dataTrigger.time{1,1};
triggerSum = sum(dataTrigger.trial{1,1},1); 
figure, plot(tm,triggerSum); hold on;
triggerSamples = find(([triggerSum>=5*10^6,0] & [0,triggerSum<=5*10^6])==1);
triggerStim01 = triggerSamples(1:3:18*3); %in block01, only response first 18 trials
scatter(tm(triggerStim01),6*10^6*ones(size(triggerStim01)),5,'filled');

%% import ECOG data and behavioral data, adapt to fieldtrip format
% load behavior data
behfile = [dRoot,'sub1_OrienWM_part02_01.mat'];
load(behfile); % offset in trial 3,4,7,17 over 20 degree; trial 19:24 no response

% ECOG data
cfg = [];
cfg.dataset = [dRoot,'direction1.edf'];
cfg.continuous = 'yes';
cfg.channel = {'all','-*DC*','-*EKG*','-*EMG*','-*BP*','-EDF*','-POL E'}; % keep channel 'POL EKG2' for later to move component EKG
dataECOG = ft_preprocessing(cfg);

fs = dataECOG.fsample;

dataECOG_block01.fsample = fs;
dataECOG_block01.label = dataECOG.label;

trialData = cell(1,3*length(triggerStim01));
trlInfoVector = nan(1,3*length(triggerStim01));
timeLine = cell(1,3*length(triggerStim01));
for t = 1:length(triggerStim01)
    trialData{1+3*(t-1)} = dataECOG.trial{1,1}(:,triggerStim01(t)-1*fs:triggerStim01(t)+3*fs);
    trialData{2+3*(t-1)} = dataECOG.trial{1,1}(:,(triggerStim01(t)-1*fs:triggerStim01(t)+3*fs)+3967); % the sample distance between item01 and item02 is 3967 
    trialData{3+3*(t-1)} = dataECOG.trial{1,1}(:,(triggerStim01(t)-1*fs:triggerStim01(t)+3*fs)+3967*2); % the sample distance between item01 and item03 is 2*3967
    trlInfoVector(1+3*(t-1)) = stim.stim01.orienChanLab(t); % orientation label (1:8 equal to 0:45:315)
    trlInfoVector(2+3*(t-1)) = stim.stim02.orienChanLab(t); % orientation label (1:8 equal to 0:45:315)
    trlInfoVector(3+3*(t-1)) = stim.stim03.orienChanLab(t); % orientation label orientation label
    timeLine{1+3*(t-1)} = -1:1/fs:3; % s
    timeLine{2+3*(t-1)} = -1:1/fs:3; % s
    timeLine{3+3*(t-1)} = -1:1/fs:3; % s
end

dataECOG_block01.trial = trialData;
dataECOG_block01.trialinfo = trlInfoVector';
dataECOG_block01.time = timeLine;

% inspect dataECOG
cfg = [];
cfg.continuous = 'no';
cfg.viewmode = 'butterfly';             % try 'butterfly','vertical'
ft_databrowser(cfg, dataECOG_block01);

% clean data, remove high variance trials
cfg = [];
cfg.method = 'summary';
dataECOG_block01_rejartif = ft_rejectvisual(cfg, dataECOG_block01);
% electrode remove: POL B30(94), POL E23(207)

% inspect data
cfg = [];
cfg.continuous = 'no';
cfg.viewmode = 'butterfly'; % also try 'butterfly','vertical'
ft_databrowser(cfg, dataECOG_block01_rejartif);

%% EKG data
cfg = [];
cfg.dataset = [dRoot,'direction1.edf'];
cfg.continuous = 'yes';
cfg.channel = {'*EKG*'};
data_EKG = ft_preprocessing(cfg);

dataEKG_block01.fsample = fs;
dataEKG_block01.label =  data_EKG.label;

trialEKGData = cell(1,3*length(triggerStim01));
for t = 1:length(triggerStim01)
    trialEKGData{1+3*(t-1)} = data_EKG.trial{1,1}(:,triggerStim01(t)-1*fs:triggerStim01(t)+3*fs);
    trialEKGData{2+3*(t-1)} = data_EKG.trial{1,1}(:,(triggerStim01(t)-1*fs:triggerStim01(t)+3*fs)+3967); % the sample distance between item01 and item02 is 3967 
    trialEKGData{3+3*(t-1)} = data_EKG.trial{1,1}(:,(triggerStim01(t)-1*fs:triggerStim01(t)+3*fs)+3967*2); % the sample distance between item01 and item03 is 2*3967
end

dataEKG_block01.trial = trialEKGData;
dataEKG_block01.trialinfo = trlInfoVector';
dataEKG_block01.time = timeLine;

% inspect data
cfg = [];
cfg.continuous = 'no';
cfg.viewmode = 'butterfly'; % also try 'butterfly','vertical'
ft_databrowser(cfg, dataEKG_block01);

%% Band pass
cfg = [];
cfg.demean = 'yes';
cfg.detrend = 'yes';
cfg.baselinewindow = [-0.2 0];
cfg.lpfilter = 'yes';
cfg.lpfreq = 200;
cfg.padding = 2;
cfg.padtype = 'data';
cfg.bsfilter = 'yes';
cfg.bsfiltord = 3;
cfg.bsfreq = [49 51; 99 101; 149 151];
dataECOG_block01_rejartif_bp = ft_preprocessing(cfg,dataECOG_block01_rejartif);
dataEKG_block01_bp = ft_preprocessing(cfg,dataEKG_block01);

% inspect data
cfg = [];
cfg.continuous = 'no';
cfg.viewmode = 'butterfly'; % also try 'butterfly','vertical'
ft_databrowser(cfg, dataECOG_block01_rejartif_bp);
ft_databrowser(cfg, dataEKG_block01_bp);

%% append ECOG and EKG data
cfg = [];
cfg.keepsampleinfo  = 'no';
data_block01_rejartif_bp_ap = ft_appenddata(cfg, dataECOG_block01_rejartif_bp, dataEKG_block01_bp);

%% downsample the data
cfg = [];
cfg.resamplefs = 500;
data_block01_rejartif_bp_ap_ds = ft_resampledata(cfg, data_block01_rejartif_bp_ap);

%% inspect and save data
cfg = [];
cfg.continuous = 'no';
cfg.viewmode = 'butterfly'; % try 'butterfly','vertical'
ft_databrowser(cfg, data_block01_rejartif_bp_ap_ds);

data_block01_preprocessed = data_block01_rejartif_bp_ap_ds; 
save('data_block01_preprocessed.mat','data_block01_preprocessed');

