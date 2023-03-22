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

% Band pass
cfg = [];
cfg.dftfilter = 'yes';
cfg.dftfreq = 50;
cfg.hpfilter='yes';
cfg.hpfreq=1;
cfg.lpfilter='yes';
cfg.lpfreq=80;

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

Part02_delay01_orientation_timefreq_item_preprocessed = data_rejartif_downsamp_comp_rejcomp_ref_item;

% load orientation channel data
subject02_part02_delay01_angle_list = [];
for i = 1:4
    loadThis = strcat(eRoot,'sub2_OrienWM_Part02_',num2str(i,'%02d'));
    load(loadThis);
    subject02_part02_delay01_angle_list = cat(1,subject02_part02_delay01_angle_list,stim.stim01.angle);
end

Part02_delay01_orientation_timefreq_item_preprocessed.trialinfo = subject02_part02_delay01_angle_list;
save('Part02_delay01_orientation_timefreq_item_preprocessed.mat','Part02_delay01_orientation_timefreq_item_preprocessed');


% Computing power spectra
cfg = [];
cfg.output = 'pow';
cfg.channel = 'all';
cfg.method  = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = 1:1:40;
cfg.toi = -0.2:0.01:1.5;
cfg.t_ftimwin = ones(size(cfg.foi)) * 0.5;
cfg.keeptrials = 'yes';

cfg.trials = find(Part02_delay01_orientation_timefreq_item_preprocessed.trialinfo == 15);
Part02_delay01_orientation_timefreq_item.orien15 = ft_freqanalysis(cfg, Part02_delay01_orientation_timefreq_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_timefreq_item_preprocessed.trialinfo == 45);
Part02_delay01_orientation_timefreq_item.orien45 = ft_freqanalysis(cfg, Part02_delay01_orientation_timefreq_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_timefreq_item_preprocessed.trialinfo == 75);
Part02_delay01_orientation_timefreq_item.orien75 = ft_freqanalysis(cfg, Part02_delay01_orientation_timefreq_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_timefreq_item_preprocessed.trialinfo == 105);
Part02_delay01_orientation_timefreq_item.orien105 = ft_freqanalysis(cfg, Part02_delay01_orientation_timefreq_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_timefreq_item_preprocessed.trialinfo == 135);
Part02_delay01_orientation_timefreq_item.orien135 = ft_freqanalysis(cfg, Part02_delay01_orientation_timefreq_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_timefreq_item_preprocessed.trialinfo == 165);
Part02_delay01_orientation_timefreq_item.orien165 = ft_freqanalysis(cfg, Part02_delay01_orientation_timefreq_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_timefreq_item_preprocessed.trialinfo == 195);
Part02_delay01_orientation_timefreq_item.orien195 = ft_freqanalysis(cfg, Part02_delay01_orientation_timefreq_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_timefreq_item_preprocessed.trialinfo == 225);
Part02_delay01_orientation_timefreq_item.orien225 = ft_freqanalysis(cfg, Part02_delay01_orientation_timefreq_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_timefreq_item_preprocessed.trialinfo == 255);
Part02_delay01_orientation_timefreq_item.orien255 = ft_freqanalysis(cfg, Part02_delay01_orientation_timefreq_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_timefreq_item_preprocessed.trialinfo == 285);
Part02_delay01_orientation_timefreq_item.orien285 = ft_freqanalysis(cfg, Part02_delay01_orientation_timefreq_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_timefreq_item_preprocessed.trialinfo == 315);
Part02_delay01_orientation_timefreq_item.orien315 = ft_freqanalysis(cfg, Part02_delay01_orientation_timefreq_item_preprocessed);

cfg.trials = find(Part02_delay01_orientation_timefreq_item_preprocessed.trialinfo == 345);
Part02_delay01_orientation_timefreq_item.orien345 = ft_freqanalysis(cfg, Part02_delay01_orientation_timefreq_item_preprocessed);

save('Part02_delay01_orientation_timefreq_item.mat','Part02_delay01_orientation_timefreq_item');

% Visualizing the power spectra
load('Part02_delay01_orientation_timefreq_item.mat')
figure;
hold on;
power15 = squeeze(mean(mean(mean(Part02_delay01_orientation_timefreq_item.orien15.powspctrm,1),2),4));
plot(Part02_delay01_orientation_timefreq_item.orien15.freq(5:40), power15(5:40), 'linewidth', 2);
power45 = squeeze(mean(mean(mean(Part02_delay01_orientation_timefreq_item.orien45.powspctrm,1),2),4));
plot(Part02_delay01_orientation_timefreq_item.orien45.freq(5:40), power45(5:40), 'linewidth', 2);
power75 = squeeze(mean(mean(mean(Part02_delay01_orientation_timefreq_item.orien75.powspctrm,1),2),4));
plot(Part02_delay01_orientation_timefreq_item.orien75.freq(5:40), power75(5:40), 'linewidth', 2);
power105 = squeeze(mean(mean(mean(Part02_delay01_orientation_timefreq_item.orien105.powspctrm,1),2),4));
plot(Part02_delay01_orientation_timefreq_item.orien105.freq(5:40), power105(5:40), 'linewidth', 2);
power135 = squeeze(mean(mean(mean(Part02_delay01_orientation_timefreq_item.orien135.powspctrm,1),2),4));
plot(Part02_delay01_orientation_timefreq_item.orien135.freq(5:40), power135(5:40), 'linewidth', 2);
power165 = squeeze(mean(mean(mean(Part02_delay01_orientation_timefreq_item.orien165.powspctrm,1),2),4));
plot(Part02_delay01_orientation_timefreq_item.orien165.freq(5:40), power165(5:40), 'linewidth', 2);
power195 = squeeze(mean(mean(mean(Part02_delay01_orientation_timefreq_item.orien195.powspctrm,1),2),4));
plot(Part02_delay01_orientation_timefreq_item.orien195.freq(5:40), power195(5:40), 'linewidth', 2);
power225 = squeeze(mean(mean(mean(Part02_delay01_orientation_timefreq_item.orien225.powspctrm,1),2),4));
plot(Part02_delay01_orientation_timefreq_item.orien225.freq(5:40), power225(5:40), 'linewidth', 2);
power255 = squeeze(mean(mean(mean(Part02_delay01_orientation_timefreq_item.orien255.powspctrm,1),2),4));
plot(Part02_delay01_orientation_timefreq_item.orien255.freq(5:40), power255(5:40), 'linewidth', 2);
power285 = squeeze(mean(mean(mean(Part02_delay01_orientation_timefreq_item.orien285.powspctrm,1),2),4));
plot(Part02_delay01_orientation_timefreq_item.orien285.freq(5:40), power285(5:40), 'linewidth', 2);
power315 = squeeze(mean(mean(mean(Part02_delay01_orientation_timefreq_item.orien315.powspctrm,1),2),4));
plot(Part02_delay01_orientation_timefreq_item.orien315.freq(5:40), power315(5:40), 'linewidth', 2);
power345 = squeeze(mean(mean(mean(Part02_delay01_orientation_timefreq_item.orien345.powspctrm,1),2),4));
plot(Part02_delay01_orientation_timefreq_item.orien345.freq(5:40), power345(5:40), 'linewidth', 2);
legend('orien15','orien45','orien75','orien105','orien135','orien165',...
    'orien195','orien225','orien255','orien285','orien315','orien345')
xlabel('Frequency (Hz)')
ylabel('Power (\mu V^2)')

% Time-frequency representations of power
cfg = [];
cfg.showlabels = 'yes';
cfg.baseline = [-0.2 0];
cfg.baselinetype = 'absolute';
cfg.layout = 'acticap-64ch-standard2.mat';
ft_multiplotTFR(cfg, Part02_delay01_orientation_timefreq_item.orien15)

% A time-frequency representation of one channel
cfg = [];
cfg.baseline = [-0.2 0];
cfg.baselinetype = 'absolute';
cfg.channel = 'Fp1';
cfg.interactive = 'no';
cfg.layout = 'acticap-64ch-standard2.mat';
ft_singleplotTFR(cfg, Part02_delay01_orientation_timefreq_item.orien15)

% A topographic representation of the time-frequency representations
cfg = [];
cfg.baseline = [-0.2 0];
cfg.baselinetype = 'absolute';
cfg.xlim = [-0.2 : 0.1 : 1.5];
cfg.ylim = [8 12];
cfg.marker = 'on';
cfg.layout = 'acticap-64ch-standard2.mat';
ft_topoplotER(cfg,Part02_delay01_orientation_timefreq_item.orien15);



