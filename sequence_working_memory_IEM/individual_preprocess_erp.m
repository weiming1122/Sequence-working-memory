clc;
clear;
close all;

restoredefaultpath                % restore default folder for matlab
maindir = pwd;                    % keep main path

cd E:\lwm\fieldtrip-master        % set up the path of fieldtrip
addpath(pwd)

ft_defaults

cd(maindir)                       % return to main
% setup directories
root = pwd; out = 'part02_delay01_analysis';
eRoot = [root(1:end-length(out)),'raw_data/'];

% Reading continuous EEG data into memory and segment
cfg = [];
cfg.dataset = [eRoot,'sub02_part02.vhdr'];

cfg.trialfun = 'ft_trialfun_general';
cfg.trialdef.eventtype = 'Stimulus';
cfg.trialdef.eventvalue = {'S242'};

cfg.trialdef.prestim = 1;
cfg.trialdef.poststim = 2;
cfg_trl_def_item = ft_definetrial(cfg);

data_raw_item = ft_preprocessing(cfg_trl_def_item);

% Baseline-correction options
cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow = [-0.2 0];

% Band pass
cfg.lpfilter='yes';
cfg.hpfilter='yes';
cfg.lpfreq=35;
cfg.hpfreq=1;

data_item = ft_preprocessing(cfg, data_raw_item);

% inspect data
cfg = [];
cfg.continuous = 'no';
cfg.viewmode = 'butterfly';             % try 'butterfly','vertical'
cfg.channel = 'EEG';
ft_databrowser(cfg, data_item);

% clean data, remove high variance trials
cfg = [];
cfg.method = 'summary';
data_rejartif_item = ft_rejectvisual(cfg, data_item);

% inspect data
cfg = [];
cfg.continuous = 'no';
cfg.viewmode = 'vertical'; % also try 'butterfly','vertical'
cfg.channel = 'EEG';
ft_databrowser(cfg, data_rejartif_item);

% downsample the data
cfg = [];
cfg.resamplefs = 200;
data_rejartif_downsamp_item = ft_resampledata(cfg, data_rejartif_item);

% perform the independent component analysis
cfg = [];
cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
data_rejartif_downsamp_comp_item = ft_componentanalysis(cfg, data_rejartif_downsamp_item);

% plot the components for visual inspection item
figure
cfg = [];
cfg.component = [1:20];       % specify the component(s) that should be plotted
cfg.layout    = 'acticap-64ch-standard2.mat'; % specify the layout file that should be used for plotting
cfg.comment   = 'no';
ft_topoplotIC(cfg, data_rejartif_downsamp_comp_item);

% For further inspection of the time course of the components item
cfg = [];
cfg.channel = [1:20]; % components to be plotted, in which we can choose which to show
cfg.layout = 'acticap-64ch-standard2.mat'; % specify the layout file that should be used for plotting
cfg.viewmode = 'component';
ft_databrowser(cfg, data_rejartif_downsamp_comp_item);

% remove the bad components and backproject the data item
cfg = [];
cfg.component = input('input the component number of removed component(s),like [n1 n2 n3]:');
data_rejartif_downsamp_comp_rejcomp_item = ft_rejectcomponent(cfg,data_rejartif_downsamp_comp_item, data_rejartif_downsamp_item);

% removing component: [1]

% inspect data
cfg = [];
cfg.continuous = 'no';
cfg.viewmode = 'butterfly'; % also try 'butterfly','vertical'
cfg.channel = 'EEG';
ft_databrowser(cfg, data_rejartif_downsamp_comp_rejcomp_item);

%reference average of all channels
cfg = [];
cfg.implicitref='FCz';
cfg.reref = 'yes';
cfg.refchannel = 'all';
cfg.refmethod = 'avg';

data_rejartif_downsamp_comp_rejcomp_ref_item = ft_preprocessing(cfg, data_rejartif_downsamp_comp_rejcomp_item);

% inspect data
cfg = [];
cfg.continuous = 'no';
cfg.viewmode = 'butterfly'; % try 'butterfly','vertical'
cfg.channel = 'EEG';
ft_databrowser(cfg, data_rejartif_downsamp_comp_rejcomp_ref_item);

Part02_delay01_orientation_ERP_item_preprocessed = data_rejartif_downsamp_comp_rejcomp_ref_item;

% load orientation channel data
subject02_part02_delay01_angle_list = [];
for i = 1:4
    loadThis = strcat(eRoot,'sub2_OrienWM_Part02_',num2str(i,'%02d'));
    load(loadThis);
    subject02_part02_delay01_angle_list = cat(1,subject02_part02_delay01_angle_list,stim.stim01.angle);
end

Part02_delay01_orientation_ERP_item_preprocessed.trialinfo = subject02_part02_delay01_angle_list;
save('Part02_delay01_orientation_ERP_item_preprocessed.mat','Part02_delay01_orientation_ERP_item_preprocessed');

% use ft_timelockanalysis to compute the ERPs
cfg = [];
cfg.keeptrials = 'yes';

cfg.trials = find(Part02_delay01_orientation_ERP_item_preprocessed.trialinfo == 15);
Part02_delay01_orientation_ERP_item.orien15 = ft_timelockanalysis(cfg, Part02_delay01_orientation_ERP_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_ERP_item_preprocessed.trialinfo == 45);
Part02_delay01_orientation_ERP_item.orien45 = ft_timelockanalysis(cfg, Part02_delay01_orientation_ERP_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_ERP_item_preprocessed.trialinfo == 75);
Part02_delay01_orientation_ERP_item.orien75 = ft_timelockanalysis(cfg, Part02_delay01_orientation_ERP_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_ERP_item_preprocessed.trialinfo == 105);
Part02_delay01_orientation_ERP_item.orien105 = ft_timelockanalysis(cfg, Part02_delay01_orientation_ERP_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_ERP_item_preprocessed.trialinfo == 135);
Part02_delay01_orientation_ERP_item.orien135 = ft_timelockanalysis(cfg, Part02_delay01_orientation_ERP_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_ERP_item_preprocessed.trialinfo == 165);
Part02_delay01_orientation_ERP_item.orien165 = ft_timelockanalysis(cfg, Part02_delay01_orientation_ERP_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_ERP_item_preprocessed.trialinfo == 195);
Part02_delay01_orientation_ERP_item.orien195 = ft_timelockanalysis(cfg, Part02_delay01_orientation_ERP_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_ERP_item_preprocessed.trialinfo == 225);
Part02_delay01_orientation_ERP_item.orien225 = ft_timelockanalysis(cfg, Part02_delay01_orientation_ERP_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_ERP_item_preprocessed.trialinfo == 255);
Part02_delay01_orientation_ERP_item.orien255 = ft_timelockanalysis(cfg, Part02_delay01_orientation_ERP_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_ERP_item_preprocessed.trialinfo == 285);
Part02_delay01_orientation_ERP_item.orien285 = ft_timelockanalysis(cfg, Part02_delay01_orientation_ERP_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_ERP_item_preprocessed.trialinfo == 315);
Part02_delay01_orientation_ERP_item.orien315 = ft_timelockanalysis(cfg, Part02_delay01_orientation_ERP_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_ERP_item_preprocessed.trialinfo == 345);
Part02_delay01_orientation_ERP_item.orien345 = ft_timelockanalysis(cfg, Part02_delay01_orientation_ERP_item_preprocessed);

save('Part02_delay01_orientation_ERP_item.mat','Part02_delay01_orientation_ERP_item');

% plots the erp for all electrodes topographically
cfg = [];
cfg.showlabels = 'yes';
cfg.layout = 'acticap-64ch-standard2.mat';
cfg.xlim = [-0.2 1.5];
cfg.baseline = [-0.2 0];
ft_multiplotER(cfg, Part02_delay01_orientation_ERP_item.orien15, Part02_delay01_orientation_ERP_item.orien45, Part02_delay01_orientation_ERP_item.orien75,...
    Part02_delay01_orientation_ERP_item.orien105, Part02_delay01_orientation_ERP_item.orien135, Part02_delay01_orientation_ERP_item.orien165,...
    Part02_delay01_orientation_ERP_item.orien195, Part02_delay01_orientation_ERP_item.orien225, Part02_delay01_orientation_ERP_item.orien255,...
    Part02_delay01_orientation_ERP_item.orien285, Part02_delay01_orientation_ERP_item.orien315, Part02_delay01_orientation_ERP_item.orien345)

% plot one sensor data use ft_singleplotER and specify the name of the channel you are interested in
cfg = [];
cfg.xlim = [-0.2 1.5];
cfg.baseline = [-0.2 0];
cfg.channel = 'Fp1';
ft_singleplotER(cfg, Part02_delay01_orientation_ERP_item.orien15, Part02_delay01_orientation_ERP_item.orien45, Part02_delay01_orientation_ERP_item.orien75,...
    Part02_delay01_orientation_ERP_item.orien105, Part02_delay01_orientation_ERP_item.orien135, Part02_delay01_orientation_ERP_item.orien165,...
    Part02_delay01_orientation_ERP_item.orien195, Part02_delay01_orientation_ERP_item.orien225, Part02_delay01_orientation_ERP_item.orien255,...
    Part02_delay01_orientation_ERP_item.orien285, Part02_delay01_orientation_ERP_item.orien315, Part02_delay01_orientation_ERP_item.orien345)

% plot the topographic distribution of the data averaged over the time interval
cfg = [];
cfg.xlim = [0 0.4];
cfg.colorbar = 'yes';
cfg.layout = 'acticap-64ch-standard2.mat';
ft_topoplotER(cfg,Part02_delay01_orientation_ERP_item.orien15);


% plot a sequence of topographic plots
cfg = [];
cfg.layout = 'acticap-64ch-standard2.mat';
cfg.xlim = [-0.2 : 0.1 : 1.5];               % Define 18 time intervals
cfg.zlim = [-0.5 1.5];                       % Set the 'color' limits.
clf;
ft_topoplotER(cfg,Part02_delay01_orientation_ERP_item.orien15);


