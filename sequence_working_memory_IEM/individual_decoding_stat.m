% compute average decoding accuracy and perform cluster-based monte carlo
% simulation analysis.
clear all;
Nsub = 1; % number of subject

Nblock = 3; % cross-validation
Nitr = 10; % iteration
Ntp = 341; % # of time points
NBins = 12; % # of location bin

tm = -200:5:1500;

%create empty matrix
AverageAccuracy = nan(Nsub,Ntp);

DecodingAccuracy = nan(Ntp,Nblock,Nitr);

%% load SVM_ECOC output files
fileLocation = pwd;
readThis =strcat(fileLocation,'/Part02_delay01_color_decoding_Results_Alphabased.mat');
load(readThis)

% prediciton from SVM-ECOC model
svmPrediction = squeeze(svmECOC.modelPredict);
tstTargets = squeeze(svmECOC.targets);
clear svmECOC

% compute decoding accuracy of each decoding trial
for block = 1:Nblock
    for itr = 1:Nitr
        for tp = 1:Ntp
            
            prediction = squeeze(svmPrediction(itr,tp,block,:)); % this is predictions from models
            TrueAnswer = squeeze(tstTargets(itr,tp,block,:)); % this is predictions from models
            Err = TrueAnswer - prediction;
            ACC = mean(Err==0);
            DecodingAccuracy(tp,block,itr) = ACC; % average decoding accuracy
            
        end
    end
end

%average across block and iterations
grandAvg = squeeze(mean(mean(DecodingAccuracy,2),3));

% perform temporal smoothing
smoothed = nan(1,Ntp);
for tAvg = 1:Ntp
    if tAvg ==1
        smoothed(tAvg) = mean(grandAvg((tAvg):(tAvg+2)));
    elseif tAvg ==2
        smoothed(tAvg) = mean(grandAvg((tAvg-1):(tAvg+2)));
    elseif tAvg == (Ntp-1)
        smoothed(tAvg) = mean(grandAvg((tAvg-2):(tAvg+1)));
    elseif tAvg == Ntp
        smoothed(tAvg) = mean(grandAvg((tAvg-2):(tAvg)));
    else
        smoothed(tAvg) = mean(grandAvg((tAvg-2):(tAvg+2)));
    end
    
end

% Save smoothe data
AverageAccuracy(Nsub,:) =smoothed; % average across iteration and block

%% plot significant clusters
figure
plot(tm,AverageAccuracy,'LineWidth',2);
xlim([-200 1500])
xlabel('Time (ms)','fontsize',14);ylabel('Decoding Accuracy','fontsize',14)
hold on
h = line(tm,1/3* ones(1,Ntp));
h.LineStyle = '--';
h.Color = [0.1,0.1,0.1];
set(gca,'linewidth',1,'fontsize',13);
title('Alpha-based decoding: color','fontsize',14)
hold off



