
function orientation_Alpha_decoding_cross_generalization

clc;
clear;
close all;

restoredefaultpath                % restore default folder for matlab
maindir = pwd;                    % keep main path

cd E:\lwm\eeglab2019_0        % set up the path of fieldtrip
addpath(genpath(pwd))

cd(maindir)                       % return to main

% parameters to set

svmECOC.nChans = 12; % # of channels
svmECOC.nBins = svmECOC.nChans; % # of stimulus bins
svmECOC.nIter = 10; % # of iterations
svmECOC.nBlocks = 3; % # of blocks for cross-validation
svmECOC.frequencies = [8 12]; % low pass filter
svmECOC.time = -200:5:1500; % time points of interest in the analysis
svmECOC.window = 2; % 1 data point per 5 ms in the preprocessed data
svmECOC.Fs = 200; % samplring rate of in the preprocessed data for filtering

ReleventChan = 1:1:64; %electrodes
svmECOC.nElectrodes = length(ReleventChan); % # of electrode included in the analysis

svmECOC.sampRate = 20; % downsampled sample rate (in ms)
svmECOC.dtime = -200:svmECOC.sampRate:1500; % downsampled time points
svmECOC.stepSize = svmECOC.sampRate/(1000/svmECOC.Fs); % number of samples the classifer jumps with each shift
svmECOC.nSamps = length(svmECOC.time); %  # of samples
svmECOC.dSamps = length(svmECOC.dtime); % # of samples at downsampled rate
% for brevity in analysis

nChans = svmECOC.nChans;
nBins = svmECOC.nBins;
nIter = svmECOC.nIter;
nBlocks = svmECOC.nBlocks;
freqs = svmECOC.frequencies;
times = svmECOC.time;
nElectrodes = svmECOC.nElectrodes;
nSamps = svmECOC.nSamps;
dtimes = svmECOC.dtime;
dSamps = svmECOC.dSamps;
stepSize = svmECOC.stepSize;
Fs = svmECOC.Fs;
window = svmECOC.window;
    
    
    % load data
    dataLocation = pwd; % set directory of data set
    loadThis = strcat(dataLocation,'/Part01_orientation_timefreq_decoding_data.mat');
    load(loadThis)
    
    % where to save decoding output
    saveLocation = pwd; % set directory for decoding results.
    
    % set up locaiton bin of each trial
    channel = Part01_orientation_timefreq_decoding_data.channel; % tip location of sample teardrop
    svmECOC.posBin = channel'; % add to fm structure so it's saved
    posBin = svmECOC.posBin;
    
    % grab EEG data
    eegs = Part01_orientation_timefreq_decoding_data.eeg(:,ReleventChan,:);
    
    % set up time points
    tois = ismember(Part01_orientation_timefreq_decoding_data.time.pre:5:Part01_orientation_timefreq_decoding_data.time.post,svmECOC.time); nTimes = length(tois);
    
    % # of trials
    svmECOC.nTrials = length(posBin); nTrials = svmECOC.nTrials;
    
    % Preallocate Matrices
    svm_predict = nan(nIter,dSamps,dSamps,nChans); % a matrix to save prediction from SVM
    tst_target = nan(nIter,dSamps,dSamps,nChans);  % a matrix to save true target values
    svmECOC.blocks = nan(nTrials,nIter);  % create svmECOC.block to save block assignments
    
    % low-pass filtering
    filtData = nan(nTrials,nElectrodes,nTimes);
    
    for c = 1:nElectrodes
        filtData(:,c,:) =abs(hilbert(eegfilt(squeeze(eegs(:,c,:)),Fs,freqs(1,1),freqs(1,2))')').^2;   %Instantaneous power
    end
    
    % trim filtered data to remove times that are not of interest (after filtering to avoid edge artifacts)
    filtData = filtData(:,:,tois);

    
    % downsample to reduced sampled rate (after filtering, so that downsampling doesn't affect filtering)
    filtData = filtData(:,:,1:stepSize:nSamps);

    
    % delay01 data
    % load data
    loadThis_delay01 = strcat(dataLocation,'/Part02_delay01_orientation_timefreq_decoding_data.mat');
    load(loadThis_delay01)
    
    % grab EEG data
    eegs_delay01 = Part02_delay01_orientation_timefreq_decoding_data.eeg(:,ReleventChan,:);
    
    channel_delay01 = Part02_delay01_orientation_timefreq_decoding_data.channel; % tip location of sample teardrop
    svmECOC.posBin_delay01 = channel_delay01'; % add to fm structure so it's saved
    posBin_delay01 = svmECOC.posBin_delay01;
    
    svmECOC.nTrials_delay01 = length(posBin_delay01); nTrials_delay01 = svmECOC.nTrials_delay01;
    % low-pass filtering
    filtData_delay01 = nan(nTrials_delay01,nElectrodes,nTimes);
    
    for c = 1:nElectrodes
        filtData_delay01(:,c,:) =abs(hilbert(eegfilt(squeeze(eegs_delay01(:,c,:)),Fs,freqs(1,1),freqs(1,2))')').^2;   %Instantaneous power
    end
    
    % trim filtered data to remove times that are not of interest (after filtering to avoid edge artifacts)
    filtData_delay01 = filtData_delay01(:,:,tois);

    
    % downsample to reduced sampled rate (after filtering, so that downsampling doesn't affect filtering)
    filtData_delay01 = filtData_delay01(:,:,1:stepSize:nSamps);

    
    % Loop through each iteration
    tic % start timing iteration loop
    for iter = 1:nIter
        
        % preallocate arrays
        blocks = nan(size(posBin));
        shuffBlocks = nan(size(posBin));
        
        % count number of trials within each position bin
        clear binCnt
        for bin = 1:nBins
            binCnt(bin) = sum(posBin == bin);
        end
        
        minCnt = min(binCnt); % # of trials for position bin with fewest trials
        nPerBin = floor(minCnt/nBlocks); % max # of trials such that the # of trials for each bin can be equated within each block
        
        % shuffle trials
        shuffInd = randperm(nTrials)'; % create shuffle index
        shuffBin = posBin(shuffInd); % shuffle trial order
        
        % take the 1st nPerBin x nBlocks trials for each position bin.
        
        for bin = 1:nBins;
            idx = find(shuffBin == bin); % get index for trials belonging to the current bin
            idx = idx(1:nPerBin*nBlocks); % drop excess trials
            x = repmat((1:nBlocks)',nPerBin,1); shuffBlocks(idx) = x; % assign randomly order trials to blocks
        end
        
        % unshuffle block assignment
        blocks(shuffInd) = shuffBlocks;
        
        % save block assignment
        svmECOC.blocks(:,iter) = blocks; % block assignment
        %        svmECOC.nTrialsPerBlock = length(blocks(blocks == 1)); % # of trials per block
        
        % Average data for each position bin across blocks
        posBins = 1:nBins;
        blockDat_filtData = nan(nBins*nBlocks,nElectrodes,dSamps);    % averaged & filtered EEG data with resampling at 50 Hz
        labels = nan(nBins*nBlocks,1);                                  % bin labels for averaged & filtered EEG data
        blockNum = nan(nBins*nBlocks,1);                                % block numbers for averaged & filtered EEG data
        bCnt = 1;
        

        for ii = 1:nBins
            for iii = 1:nBlocks
                blockDat_filtData(bCnt,:,:) = squeeze(mean(filtData(posBin==posBins(ii) & blocks==iii,:,:),1));
                labels(bCnt) = ii;
                blockNum(bCnt) = iii;
                bCnt = bCnt+1;
            end
        end
        
        % preallocate arrays
        blocks_delay01 = nan(size(posBin_delay01));
        shuffBlocks_delay01 = nan(size(posBin_delay01));
        
        clear binCnt_delay01
        for bin = 1:nBins
            binCnt_delay01(bin) = sum(posBin_delay01 == bin);
        end
        minCnt_delay01 = min(binCnt_delay01);
        if minCnt_delay01 >= nPerBin
            minCnt_delay01 = nPerBin;
        end
        nPerBin_delay01 = minCnt_delay01;
        % take the 1st nPerBin_delay01 trials for each position bin.
        % shuffle trials
        shuffInd_delay01 = randperm(nTrials_delay01)'; % create shuffle index
        shuffBin_delay01 = posBin_delay01(shuffInd_delay01); % shuffle trial order
        
        for bin = 1:nBins;
            idx = find(shuffBin_delay01 == bin); % get index for trials belonging to the current bin
            idx = idx(1:nPerBin_delay01); % drop excess trials
            shuffBlocks_delay01(idx) = nBlocks+1;
        end
        
        % unshuffle block assignment
        blocks_delay01(shuffInd_delay01) = shuffBlocks_delay01;
        
        % save block assignment
        svmECOC.blocks_delay01(:,iter) = blocks_delay01; % block assignment
        
        % Average data for each position bin across blocks
        blockDat_filtData_delay01 = nan(nBins,nElectrodes,dSamps);    % averaged & filtered EEG data with resampling at 50 Hz
        labels_delay01 = nan(nBins,1);                                  % bin labels for averaged & filtered EEG data
        blockNum_delay01 = nan(nBins,1);                                % block numbers for averaged & filtered EEG data
        bCnt_delay01 = 1;
        
        for ii = 1:nBins
            for iii = 4
                blockDat_filtData_delay01(bCnt_delay01,:,:) = squeeze(mean(filtData_delay01(posBin_delay01==posBins(ii) & blocks_delay01==iii,:,:),1));
                labels_delay01(bCnt_delay01) = ii;
                blockNum_delay01(bCnt_delay01) = iii;
                bCnt_delay01 = bCnt_delay01+1;
            end
        end
        
        % Do SVM_ECOC at each time point
        for t = 1:dSamps
%             % grab data for timepoint t
%             toi = ismember(times,times(t)-window/2:times(t)+window/2);
%             % average across time window of interest
%             dataAtTimeT = squeeze(mean(blockDat_filtData(:,:,toi),3));
%             dataAtTimeT_delay01 = squeeze(mean(blockDat_filtData_delay01(:,:,toi),3));
%             toi = ismember(times,times(t)-window/2:times(t)+window/2);
%             dataAtTimeT = squeeze(mean(blockDat_filtData(:,:,[find(times==400):find(times==600)]),3));
%             dataAtTimeT_delay01 = squeeze(mean(blockDat_filtData_delay01(:,:,toi),3));
            
            toi = ismember(dtimes,dtimes(t)-window/2:dtimes(t)+window/2);
            dataAtTimeT = squeeze(mean(blockDat_filtData(:,:,toi),3));
            for ii = 1: dSamps
            dataAtTimeT_delay01 = squeeze(mean(blockDat_filtData_delay01(:,:,ii),3));
            
            % Do SVM_ECOC for each block
            trnl = labels(blockNum~=4); % training labels
            tstl = labels_delay01(blockNum_delay01==4); % test labels
            trnD = dataAtTimeT(blockNum~=4,:);    % training data
            tstD = dataAtTimeT_delay01(blockNum_delay01==4,:);    % test data
            
            % SVM_ECOC
            mdl = fitcecoc(trnD,trnl, 'Coding','onevsall','Learners','SVM' );   %train support vector mahcine
            LabelPredicted = predict(mdl, tstD);       % predict classes for new data
            svm_predict(iter,t,ii,:) = LabelPredicted;  % save predicted labels
            tst_target(iter,t,ii,:) = tstl;             % save true target labels
            end
        end % end of time points
    end % end of iteration
    toc % stop timing the iteration loop
    
    OutputfName = strcat(saveLocation,'/Part02_delay01_orientation_decoding_Results_Alphabased_generalization.mat');
    
    svmECOC.targets = tst_target;
    svmECOC.modelPredict = svm_predict;
    svmECOC.nBlocks = nBlocks;
    save(OutputfName,'svmECOC','-v7.3');


