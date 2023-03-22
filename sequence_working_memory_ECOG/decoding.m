function sub01_decoding_orientation_Alpha(subs)
clc;
clear;
close all;

delete(gcp)
parpool

maindir = pwd;                    % keep main path

cd E:\lwm\eeglab2019_0        % set up the path of fieldtrip
addpath(genpath(pwd))

cd(maindir)                       % return to main

% participants ID
if nargin ==0
    subs = [1];
end
nSubs = length(subs);

% parameters to set
svmECOC.nChans = 8; % # of channels
svmECOC.nBins = svmECOC.nChans; % # of stimulus bins
svmECOC.nIter = 10; % # of iterations
svmECOC.nBlocks = 3; % # of blocks for cross-validation
svmECOC.frequencies = [8 12]; % low pass filter
svmECOC.time = -0.200:0.002:2.000; % time points of interest in the analysis
svmECOC.alltime = -1.000:0.002:3.000; % time points of interest in the analysis
svmECOC.window = 2; % 1 data point per 5 ms in the preprocessed data
svmECOC.Fs = 200; % samplring rate of in the preprocessed data for filtering

ReleventChan = 1:1:206; %electrodes
svmECOC.nElectrodes = length(ReleventChan); % # of electrode included in the analysis

% for brevity in analysis
nChans = svmECOC.nChans;
nBins = svmECOC.nBins;
nIter = svmECOC.nIter;
nBlocks = svmECOC.nBlocks;
freqs = svmECOC.frequencies;
times = svmECOC.time;
alltimes = svmECOC.alltime;
nElectrodes = svmECOC.nElectrodes;
nSamps = length(svmECOC.time);
Fs = svmECOC.Fs;

%% Loop through participants
for s = 1:nSubs
    sn = subs(s);
    fprintf('Subject:\t%d\n',sn)
    % load data
    currentSub = num2str(sn,'%02d');
    dataLocation = pwd; % set directory of data set
    loadThis = strcat(dataLocation,'/subject',currentSub,'_DecodingData.mat');
    load(loadThis);
    
    % where to save decoding output
    saveLocation = pwd; % set directory for decoding results.
    
    % set up locaiton bin of each trial
    channel = DecodingData.trialinfo; % tip location of sample teardrop
    svmECOC.posBin = channel'; % add to fm structure so it's saved
    posBin = svmECOC.posBin;
    
    % grab EEG data
    for i= 1:length(posBin)
        eegs(i,:,:) = DecodingData.trial{1,i}(ReleventChan,:);
    end
    
    % set up time points
    tois = ismember(round(1000*alltimes),round(1000*times)); nTimes = length(tois);
    
    % # of trials
    svmECOC.nTrials = length(posBin); nTrials = svmECOC.nTrials;
    
    % Preallocate Matrices
    svm_predict = nan(nIter,nSamps,nBlocks,nChans); % a matrix to save prediction from SVM
    tst_target = nan(nIter,nSamps,nBlocks,nChans);  % a matrix to save true target values
    svmECOC.blocks = nan(nTrials,nIter);  % create svmECOC.block to save block assignments
    
    % low-pass filtering
    filtData = nan(nTrials,nElectrodes,nTimes);
    
    for c = 1:nElectrodes
        filtData(:,c,:) = eegfilt(squeeze(eegs(:,c,:)),Fs,freqs(1,1),freqs(1,2)); % low pass filter
    end
    
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
        for bin = 1:nBins
            idx = find(shuffBin == bin); % get index for trials belonging to the current bin
            idx = idx(1:nPerBin*nBlocks); % drop excess trials
            x = repmat((1:nBlocks)',nPerBin,1); shuffBlocks(idx) = x; % assign randomly order trials to blocks
        end
        
        % unshuffle block assignment
        blocks(shuffInd) = shuffBlocks;
        
        % save block assignment
        svmECOC.blocks(:,iter) = blocks; % block assignment
        svmECOC.nTrialsPerBlock = length(blocks(blocks == 1)); % # of trials per block
        
        % Average data for each position bin across blocks
        posBins = 1:nBins;
        blockDat_filtData = nan(nBins*nBlocks,nElectrodes,nSamps);    % averaged & filtered EEG data with resampling at 50 Hz
        labels = nan(nBins*nBlocks,1);                                  % bin labels for averaged & filtered EEG data
        blockNum = nan(nBins*nBlocks,1);                                % block numbers for averaged & filtered EEG data
        bCnt = 1;
        
        for ii = 1:nBins
            for iii = 1:nBlocks
                blockDat_filtData(bCnt,:,:) = squeeze(mean(filtData(posBin==posBins(ii) & blocks==iii,:,tois),1));
                labels(bCnt) = ii;
                blockNum(bCnt) = iii;
                bCnt = bCnt+1;
            end
        end
        
        % Do SVM_ECOC at each time point
        parfor t = 1:nSamps
            % grab data for timepoint t
            toi = ismember(times,times(t));
            % average across time window of interest
            dataAtTimeT = squeeze(mean(blockDat_filtData(:,:,toi),3));
            
            % Do SVM_ECOC for each block
            for i=1:nBlocks % loop through blocks, holding each out as the test set
                trnl = labels(blockNum~=i); % training labels
                tstl = labels(blockNum==i); % test labels
                trnD = dataAtTimeT(blockNum~=i,:);    % training data
                tstD = dataAtTimeT(blockNum==i,:);    % test data
                
                % SVM_ECOC
                mdl = fitcecoc(trnD,trnl, 'Coding','onevsall','Learners','SVM' );   %train support vector mahcine
                LabelPredicted = predict(mdl, tstD);       % predict classes for new data
                svm_predict(iter,t,i,:) = LabelPredicted;  % save predicted labels
                tst_target(iter,t,i,:) = tstl;             % save true target labels
                
            end % end of block
        end % end of time points
    end % end of iteration
    toc % stop timing the iteration loop
    
    OutputfName = strcat(saveLocation,'/Orientation_Results_Alphabased_',currentSub,'.mat');
    svmECOC.targets = tst_target;
    svmECOC.modelPredict = svm_predict;
    svmECOC.nBlocks = nBlocks;
    save(OutputfName,'svmECOC','-v7.3');
    
end % end of subject loop
