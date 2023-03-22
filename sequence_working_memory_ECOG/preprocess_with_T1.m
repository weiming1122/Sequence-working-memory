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

file_mr = [dRoot,'ycb_mr_t1.nii'];
file_ct = [dRoot,'HEAD_1_0_H31S_0002_03_Head_1mm_20191125203407_2.nii'];
file_ecog = [dRoot,'direction1.edf'];

%% Anatomical workflow
subjID = 'sub01';

%% Preprocessing of the anatomical MRI

% Import the anatomical MRI
mri = ft_read_mri(file_mr);  % or a single file of a DICOM series

% Determine the native orientation of the anatomical MRI¡¯s left-right axis
ft_determine_coordsys(mri); % x goes from left to right (positive values at left side of the head)

% Align the anatomical MRI to the ACPC coordinate system
% Selecting, what is the right hemisphere (R), anterior (A), posterior comisure
% (P) and top part of the brain (Z). Quit (Q)
cfg = [];
cfg.method = 'interactive'; 
cfg.coordsys = 'acpc'; % anterior & posterior commissure, positive midline, and right side
mri_acpc = ft_volumerealign(cfg, mri);

% Write the preprocessed anatomical MRI out to file, for later processing
cfg = [];
cfg.filename = [subjID '_MR_acpc'];
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, mri_acpc);

%% Running freesurfer parcellation
% Execute FreeSurfer¡¯s recon-all functionality from the Linux (~10 hr)
% then copy results data into raw_data folder
fshome = 'D:\linuxToolbox\freesurfer';
subdir = 'D:\linuxToolbox\sub01';
mrfile  = 'D:\linuxToolbox\ycb_mr_t1.nii';
'C:\Windows\System32\bash.exe';
system(['export FREESURFER_HOME=' fshome '; ' ...
'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
'mri_convert -c -oc 0 0 0 ' mrfile ' ' [subdir '/tmp.nii'] '; ' ...
'recon-all -i ' [subdir '/tmp.nii'] ' -s ' 'freesurfer' ' -sd ' subdir ' -all'])

% visualize the pial surface of the left hemisphere
pial_lh = ft_read_headshape([dRoot,'sub01_freesurfer/surf/lh.pial']);
pial_lh.coordsys = 'acpc';
ft_plot_mesh(pial_lh);
lighting gouraud;
camlight;
% visualize the pial surface of the right hemisphere
pial_rh = ft_read_headshape([dRoot,'sub01_freesurfer/surf/rh.pial']);
pial_rh.coordsys = 'acpc';
ft_plot_mesh(pial_rh);
lighting gouraud;
camlight;

% Import the FreeSurfer-processed MRI
fsmri_acpc = ft_read_mri([dRoot,'sub01_freesurfer/mri/T1.nii']); %
fsmri_acpc.coordsys = 'acpc';

%% Preprocessing of the anatomical CT and co-register to the MRI

% Import the anatomical CT
ct = ft_read_mri(file_ct); % or a single file of a DICOM series

% Determine the native orientation of the anatomical ct's left-right axis
ft_determine_coordsys(ct);

% Align the anatomical CT to the CTF head surface coordinate system
cfg = [];
cfg.method = 'interactive'; % nasion, left & right ear, and positive midline
cfg.coordsys = 'ctf';
ct_ctf = ft_volumerealign(cfg, ct);

% Automatically convert the CT¡¯s coordinate system into an approximation of the ACPC coordinate system
ct_acpc = ft_convert_coordsys(ct_ctf, 'acpc',  [], 0); % convert ct to spm/acpac coordsys
% if you dont have spm installed it seems that you need to use method 0

% Fuse the CT with the MRI
% link the electrode locations in the anatomical CT to their corresponding locations in the anatomical MRI
cfg = [];
cfg.method = 'spm';
cfg.spmversion = 'spm12';
cfg.coordsys = 'acpc';
cfg.viewresult = 'yes'; 
ct_acpc_f = ft_volumerealign(cfg, ct_acpc, fsmri_acpc); % co-register the CT to the MRI

% Write the MRI-fused anatomical CT out to file, as a backup
cfg = [];
cfg.filename = [subjID '_CT_acpc_f'];
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, ct_acpc_f);

% check fuse result
ft_determine_coordsys(ct_acpc_f);

%% Electrode placement
hdr = ft_read_header(file_ecog);

% Localize the electrodes in the post-implant CT
cfg = [];
cfg.channel = hdr.label;
elec_acpc_f = ft_electrodeplacement(cfg, ct_acpc_f, fsmri_acpc);

% check the chanpos field
elec_acpc_f

% Visualize the MRI along with the electrodes and their labels
ft_plot_ortho(fsmri_acpc.anatomy, 'transform', fsmri_acpc.transform, 'style', 'intersect');
ft_plot_sens(elec_acpc_f, 'label', 'on', 'fontcolor', 'w');

save([subjID '_elec_acpc_f.mat'], 'elec_acpc_f');

% Brain shift compensation (optional for cortical grids and strips)
cfg = [];
cfg.method = 'cortexhull';
cfg.headshape = [dRoot,'sub01_freesurfer/surf/lh.pial'];
cfg.fshome = 'D:\linuxToolbox\freesurfer'; 
hull_lh = ft_prepare_mesh(cfg);

% Visualize the cortex and electrodes together
ft_plot_mesh(pial_lh);
ft_plot_mesh(pial_rh);
ft_plot_sens(elec_acpc_f);
view([-55 10]);
material dull;
lighting gouraud;
camlight;

% Volume-based registration (optional)
cfg = [];
cfg.nonlinear = 'yes';
cfg.spmversion = 'spm12';
cfg.spmmethod = 'new';
fsmri_mni = ft_volumenormalise(cfg, fsmri_acpc);

elec_mni_frv = elec_acpc_f;
elec_mni_frv.elecpos = ft_warp_apply(fsmri_mni.params, elec_acpc_f.elecpos, 'individual2sn');
elec_mni_frv.chanpos = ft_warp_apply(fsmri_mni.params, elec_acpc_f.chanpos, 'individual2sn');
elec_mni_frv.coordsys = 'mni';

% Visualize the cortical mesh extracted from the standard MNI brain 
% along with the spatially normalized electrodes
[ftver, ftpath] = ft_version;
load([ftpath filesep 'template/anatomy/surface_pial_left.mat']);
ft_plot_mesh(mesh);
ft_plot_sens(elec_mni_frv);
view([-90 20]);
material dull;
lighting gouraud;
camlight;

% Save the normalized electrode information to file.
save([subjID '_elec_mni_frv.mat'], 'elec_mni_frv');

