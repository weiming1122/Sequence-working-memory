clc;
clear;
close all;

%% setup directories
root = pwd; out = 'data_analysis';
dRoot = [root(1:end-length(out)),'raw_data/'];

%% load preprocessed ECOG data
load('data_block01_preprocessed.mat');
load('data_block02_preprocessed.mat');
load('data_block03_preprocessed.mat');
load('data_block04_preprocessed.mat');

%% seletdata to have same channel label for 4 blocks
cfg = [];
cfg.channel = {'all','-POL B30','-POL E23','-POL A23','-POL EKG1'};
data_block01_preprocessed_equal = ft_selectdata(cfg, data_block01_preprocessed);
data_block02_preprocessed_equal = ft_selectdata(cfg, data_block02_preprocessed);
data_block03_preprocessed_equal = ft_selectdata(cfg, data_block03_preprocessed);
data_block04_preprocessed_equal = ft_selectdata(cfg, data_block04_preprocessed);

%% append data
cfg = [];
cfg.keepsampleinfo  = 'no';
data_preprocessed_equal = ft_appenddata(cfg, data_block01_preprocessed_equal,data_block02_preprocessed_equal, ...
    data_block03_preprocessed_equal, data_block04_preprocessed_equal);

% inspect data
cfg = [];
cfg.continuous = 'no';
cfg.viewmode = 'butterfly'; % try 'butterfly','vertical'
ft_databrowser(cfg, data_preprocessed_equal);

%% perform the independent component analysis
cfg = [];
cfg.method = 'runica';
data_preprocessed_equal_comp = ft_componentanalysis(cfg, data_preprocessed_equal);

save('subject01_data_preprocessed_equal_comp.mat','data_preprocessed_equal_comp','-v7.3');

% prepare for lay
pial_lh = ft_read_headshape([dRoot,'sub01_freesurfer/surf/lh.pial']);
pial_lh.coordsys = 'acpc';
cfg = [];
cfg.headshape = pial_lh;
cfg.projection = 'orthographic';
cfg.channel = {'all'};
cfg.viewpoint  = 'left';
cfg.elec = 'sub01_elec_acpc_f.mat';
lay = ft_prepare_layout(cfg);

% plot the components for visual inspection item
figure
cfg = [];
cfg.component = [1:40];       % specify the component(s) that should be plotted
cfg.layout = lay; % specify the layout file that should be used for plotting
cfg.comment = 'no';
ft_topoplotIC(cfg, data_preprocessed_equal_comp);

% For further inspection of the time course of the components
cfg = [];
cfg.component = [1:40];       % components to be plotted, in which we can choose which to show
cfg.layout = lay; % specify the layout file that should be used for plotting
cfg.viewmode = 'component';
ft_databrowser(cfg, data_preprocessed_equal_comp);

% remove the bad components and backproject the data
cfg = [];
cfg.component = input('input the component number of removed component(s),like [n1 n2 n3]:');
data_preprocessed_equal_comp_rejcomp = ft_rejectcomponent(cfg,...
    data_preprocessed_equal_comp, data_preprocessed_equal);
% removing component: [11 18 20 39 75 83 85]

%%
% % Correlate EKG with ICs
% ekg_ic_corr = zeros([numel(data_preprocessed_equal_comp.topolabel), numel(data_preprocessed_equal_comp.trial)]);
% for ic_ix = 1:numel(data_preprocessed_equal_comp.topolabel)
%     for t_ix = 1:numel(dataECOG_block01_rejartif_bp_ds_comp.trial)
%         temp = corrcoef(dataEKG_block01_bp_ds.trial{t_ix}(ekg_ix,:), dataECOG_block01_rejartif_bp_ds_comp.trial{t_ix}(ic_ix,:));
%         ekg_ic_corr(ic_ix, t_ix) = temp(1,2);
%     end
% end
% avg_ekg_ic_corr = mean(ekg_ic_corr,2);

DecodingData = data_preprocessed_equal_comp_rejcomp;

% inspect data
cfg = [];
cfg.continuous = 'no';
cfg.viewmode = 'butterfly';             % try 'butterfly','vertical'
ft_databrowser(cfg, DecodingData);

save('subject01_DecodingData.mat','DecodingData','-v7.3');

%% Time-frequency Analysis
clc;
clear;
close all;
load('subject01_DecodingData.mat');

%% Computing power spectra
cfg = [];
cfg.output = 'pow';
cfg.channel = 'all';
cfg.method  = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = 5:5:200;
cfg.toi = -0.2:0.01:2.0;
cfg.t_ftimwin = ones(size(cfg.foi)) * 0.2;
cfg.keeptrials = 'no';
DecodingData_timefreq = ft_freqanalysis(cfg, DecodingData);

cfg = [];
cfg.baseline = [-0.2 0];
cfg.baselinetype = 'relchange'; % 'absolute', 'relative', 'relchange', 'normchange', 'db', 'vssum' or 'zscore' 
DecodingData_timefreq_blc = ft_freqbaseline(cfg, DecodingData_timefreq);

cfg = [];
cfg.baseline = [-0.2 0];
cfg.baselinetype = 'db'; % 'absolute', 'relative', 'relchange', 'normchange', 'db', 'vssum' or 'zscore' 
DecodingData_timefreq_db = ft_freqbaseline(cfg, DecodingData_timefreq);

%% plot power spectra
figure;
hold on;
power = squeeze(mean(mean(DecodingData_timefreq_blc.powspctrm,1),3));
plot(DecodingData_timefreq_blc.freq(1:40), power(1:40), 'linewidth', 2);
xlabel('Frequency (Hz)','fontsize',14); ylabel('Relative Change','fontsize',14);
set(gca,'linewidth',1,'fontsize',13);
title('power spectra (relative)','fontsize',14);
hold off

figure;
hold on;
power = squeeze(mean(mean(DecodingData_timefreq_db.powspctrm,1),3));
plot(DecodingData_timefreq_db.freq(1:40), power(1:40), 'linewidth', 2);
xlabel('Frequency (Hz)','fontsize',14); ylabel('Relative Change','fontsize',14);
set(gca,'linewidth',1,'fontsize',13);
title('power spectra (relative)','fontsize',14);
hold off

figure;
hold on;
power = squeeze(mean(mean(DecodingData_timefreq.powspctrm,1),3));
plot(DecodingData_timefreq_blc.freq(1:40), power(1:40), 'linewidth', 2);
xlabel('Frequency (Hz)','fontsize',14); ylabel('Power (\mu V^2)','fontsize',14);
set(gca,'linewidth',1,'fontsize',13);
title('power spectra (absolute)','fontsize',14);
hold off

%% Time-frequency representations of power blc
%% Plotting single average across electrodes of Central (1)
cfg = [];
cfg.channel = {'*A*'};
cfg.colormap = parula;
cfg.xlim = [-0.2 2.0];
cfg.ylim = [5 200];
figure(1);
ft_singleplotTFR(cfg,DecodingData_timefreq_blc);
set(gca,'Fontsize',14);
title('Mean over Central electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('frequency (Hz)');
c = colorbar;
c.LineWidth = 1;
c.FontSize  = 12;
title(c,'Rel. change');
saveas(figure(1),'time_frequency_blc_central','png')

%% Plotting single average across electrodes of Frontal (2)
cfg = [];
cfg.channel = {'*B*'};
cfg.colormap = parula;
cfg.xlim = [-0.2 2.0];
cfg.ylim = [5 200];
figure(2);
ft_singleplotTFR(cfg,DecodingData_timefreq_blc);
set(gca,'Fontsize',14);
title('Mean over Frontal electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('frequency (Hz)');
c = colorbar;
c.LineWidth = 1;
c.FontSize  = 12;
title(c,'Rel. change');
saveas(figure(2),'time_frequency_blc_frontal','png')

%% Plotting single average across electrodes of Temporal (3)
cfg = [];
cfg.channel = {'*C*'};
cfg.colormap = parula;
cfg.xlim = [-0.2 2.0];
cfg.ylim = [5 200];
figure(3);
ft_singleplotTFR(cfg,DecodingData_timefreq_blc);
set(gca,'Fontsize',14);
title('Mean over Temporal electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('frequency (Hz)');
c = colorbar;
c.LineWidth = 1;
c.FontSize  = 12;
title(c,'Rel. change');
saveas(figure(3),'time_frequency_blc_temporal','png')

%% Plotting single average across electrodes of Wernicke (4)
cfg = [];
cfg.channel = {'*D*'};
cfg.colormap = parula;
cfg.xlim = [-0.2 2.0];
cfg.ylim = [5 200];
figure(4);
ft_singleplotTFR(cfg,DecodingData_timefreq_blc);
set(gca,'Fontsize',14);
title('Mean over Wernicke electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('frequency (Hz)');
c = colorbar;
c.LineWidth = 1;
c.FontSize  = 12;
title(c,'Rel. change');
saveas(figure(4),'time_frequency_blc_wernicke','png')

%% Plotting single average across electrodes of Broca (5)
cfg = [];
cfg.channel = {'*E*'};
cfg.colormap = parula;
cfg.xlim = [-0.2 2.0];
cfg.ylim = [5 200];
figure(5);
ft_singleplotTFR(cfg,DecodingData_timefreq_blc);
set(gca,'Fontsize',14);
title('Mean over Broca electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('frequency (Hz)');
c = colorbar;
c.LineWidth = 1;
c.FontSize  = 12;
title(c,'Rel. change');
saveas(figure(5),'time_frequency_blc_broca','png')

% Multi plot
cfg = [];
cfg.layout = lay;
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
cfg.showlabels = 'yes';
cfg.fontsize = 10;
cfg.comment = sprintf('\n');
cfg.xlim = [-0.2 2.0];
cfg.ylim = [5 200];
cfg.colorbar = 'yes'; % or 'southoutside'
cfg.colormap = parula; % 'parula' or 'jet'
figure
ft_multiplotTFR(cfg, DecodingData_timefreq_blc);
title('time frequency (relchange)');

%% Time-frequency representations of power db
%% Plotting single average across electrodes of Central (6)
cfg = [];
cfg.channel = {'*A*'};
cfg.colormap = parula;
cfg.xlim = [-0.2 2.0];
cfg.ylim = [5 200];
figure(6);
ft_singleplotTFR(cfg,DecodingData_timefreq_db);
set(gca,'Fontsize',14);
title('Mean over Central electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('frequency (Hz)');
c = colorbar;
c.LineWidth = 1;
c.FontSize  = 12;
title(c,'db');
saveas(figure(6),'time_frequency_db_central','png')

%% Plotting single average across electrodes of Frontal (7)
cfg = [];
cfg.channel = {'*B*'};
cfg.colormap = parula;
cfg.xlim = [-0.2 2.0];
cfg.ylim = [5 200];
figure(7);
ft_singleplotTFR(cfg,DecodingData_timefreq_db);
set(gca,'Fontsize',14);
title('Mean over Frontal electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('frequency (Hz)');
c = colorbar;
c.LineWidth = 1;
c.FontSize  = 12;
title(c,'db');
saveas(figure(7),'time_frequency_db_frontal','png')

%% Plotting single average across electrodes of Temporal (8)
cfg = [];
cfg.channel = {'*C*'};
cfg.colormap = parula;
cfg.xlim = [-0.2 2.0];
cfg.ylim = [5 200];
figure(8);
ft_singleplotTFR(cfg,DecodingData_timefreq_db);
set(gca,'Fontsize',14);
title('Mean over Temporal electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('frequency (Hz)');
c = colorbar;
c.LineWidth = 1;
c.FontSize  = 12;
title(c,'db');
saveas(figure(8),'time_frequency_db_temporal','png')

%% Plotting single average across electrodes of Wernicke (9)
cfg = [];
cfg.channel = {'*D*'};
cfg.colormap = parula;
cfg.xlim = [-0.2 2.0];
cfg.ylim = [5 200];
figure(9);
ft_singleplotTFR(cfg,DecodingData_timefreq_db);
set(gca,'Fontsize',14);
title('Mean over Wernicke electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('frequency (Hz)');
c = colorbar;
c.LineWidth = 1;
c.FontSize  = 12;
title(c,'db');
saveas(figure(9),'time_frequency_db_wernicke','png')

%% Plotting single average across electrodes of Broca (10)
cfg = [];
cfg.channel = {'*E*'};
cfg.colormap = parula;
cfg.xlim = [-0.2 2.0];
cfg.ylim = [5 200];
figure(10);
ft_singleplotTFR(cfg,DecodingData_timefreq_db);
set(gca,'Fontsize',14);
title('Mean over Broca electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('frequency (Hz)');
c = colorbar;
c.LineWidth = 1;
c.FontSize  = 12;
title(c,'db');
saveas(figure(10),'time_frequency_db_broca','png')

% Multi plot
cfg = [];
cfg.layout = lay;
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
cfg.showlabels = 'yes';
cfg.fontsize = 10;
cfg.comment = sprintf('\n');
cfg.xlim = [-0.2 2.0];
cfg.ylim = [5 200];
cfg.colorbar = 'yes'; % or 'southoutside'
cfg.colormap = parula; % 'parula' or 'jet'
figure
ft_multiplotTFR(cfg, DecodingData_timefreq_db);
title('time frequency (db)');

%% ERP Analysis
clc;
clear;
close all;
load('subject01_DecodingData.mat');

%% Calculate ERPs over Central electrodes (1)
cfg = [];
cfg.channel = {'*A*'};
% cfg.trials = [1:3:252];
erpCentral = ft_timelockanalysis(cfg,DecodingData);

% Plotting mean ERPs over Central electrodes
cfg = [];
cfg.channel = {'*A*'};
cfg.linewidth = 2;
cfg.graphcolor = 'br';
figure(1);
ft_singleplotER(cfg,erpCentral);
legend({'ERP'})
% Matlab commands to make the figure intelligible
set(gca,'Fontsize',14);
title('Mean over central electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('amplitude (\muV)');
h = line(-1:0.01:3,0*ones(1,length([-1:0.01:3])));
h.LineStyle = '--';
h.Color = [0.1,0.1,0.1];
saveas(figure(1),'erp_central','png')

%% Calculate ERPs over Frontal electrodes (2)
cfg = [];
cfg.channel = {'*B*'};
% cfg.trials = [1:3:252];
erpFrontal = ft_timelockanalysis(cfg,DecodingData);

% Plotting mean ERPs over Frontal electrodes
cfg = [];
cfg.channel = {'*B*'};
cfg.linewidth = 2;
cfg.graphcolor = 'br';
figure(2);
ft_singleplotER(cfg,erpFrontal);
legend({'ERP'})
% Matlab commands to make the figure intelligible
set(gca,'Fontsize',14);
title('Mean over frontal electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('amplitude (\muV)');
h = line(-1:0.01:3,0*ones(1,length([-1:0.01:3])));
h.LineStyle = '--';
h.Color = [0.1,0.1,0.1];
saveas(figure(2),'erp_frontal','png')

%% Calculate ERPs over Temporal electrodes (3)
cfg = [];
cfg.channel = {'*C*'};
% cfg.trials = [1:3:252];
erpTemporal = ft_timelockanalysis(cfg,DecodingData);

% Plotting mean ERPs over Temporal electrodes
cfg = [];
cfg.channel = {'*C*'};
cfg.linewidth = 2;
cfg.graphcolor = 'br';
figure(3);
ft_singleplotER(cfg,erpTemporal);
legend({'ERP'})
% Matlab commands to make the figure intelligible
set(gca,'Fontsize',14);
title('Mean over temporal electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('amplitude (\muV)');
h = line(-1:0.01:3,0*ones(1,length([-1:0.01:3])));
h.LineStyle = '--';
h.Color = [0.1,0.1,0.1];
saveas(figure(3),'erp_temporal','png')

%% Calculate ERPs over Wernicke electrodes (4)
cfg = [];
cfg.channel = {'*D*'};
% cfg.trials = [1:3:252];
erpWernicke = ft_timelockanalysis(cfg,DecodingData);

% Plotting mean ERPs over Wernicke electrodes
cfg = [];
cfg.channel = {'*D*'};
cfg.linewidth = 2;
cfg.graphcolor = 'br';
figure(4);
ft_singleplotER(cfg,erpWernicke);
legend({'ERP'})
% Matlab commands to make the figure intelligible
set(gca,'Fontsize',14);
title('Mean over wernicke electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('amplitude (\muV)');
h = line(-1:0.01:3,0*ones(1,length([-1:0.01:3])));
h.LineStyle = '--';
h.Color = [0.1,0.1,0.1];
saveas(figure(4),'erp_wernicke','png')

%% Calculate ERPs over Broca electrodes (5)
cfg = [];
cfg.channel = {'*E*'};
% cfg.trials = [1:3:252];
erpBroca = ft_timelockanalysis(cfg,DecodingData);

% Plotting mean ERPs over Broca electrodes
cfg = [];
cfg.channel = {'*E*'};
cfg.linewidth = 2;
cfg.graphcolor = 'br';
figure(5);
ft_singleplotER(cfg,erpBroca);
legend({'ERP'})
% Matlab commands to make the figure intelligible
set(gca,'Fontsize',14);
title('Mean over broca electrodes');
set(gca,'box','on');
xlabel('time (s)');
ylabel('amplitude (\muV)');
h = line(-1:0.01:3,0*ones(1,length([-1:0.01:3])));
h.LineStyle = '--';
h.Color = [0.1,0.1,0.1];
saveas(figure(5),'erp_broca','png')
