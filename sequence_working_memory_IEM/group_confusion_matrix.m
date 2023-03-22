clc;
clear;
close all;

maindir = pwd;                    % keep main path
cd E:\lwm\plotConfMat-master
addpath(pwd)
cd(maindir)                       % return to main

subList = [2:9];
Nsub = length(subList);                % return to main

% Loop through participants
for s = 1:Nsub
    sn = subList(s);
    fprintf('Subject:\t%d\n',sn)
    
    fileLocation = pwd;
    readThis =strcat(fileLocation,'/Color_Results_Alphabased_',num2str(sn,'%02d'),'.mat')
    load(readThis)
    
    for i = 1:length(svmECOC.time)
        a = squeeze(svmECOC.targets(:,i,:,:));
        b = squeeze(svmECOC.modelPredict(:,i,:,:));
        raw= reshape(a, [3*10*3,1]);
        predict = reshape(b, [3*10*3,1]);
        [matrix,order] = confusionmat(raw,predict);
        confusion(:,:,i,s) = matrix;
    end
end

for nTime = 61
    timepoint_start = svmECOC.time(nTime);
    timepoint_end = svmECOC.time(nTime+60);
    confusion_averagetime = squeeze(mean(mean(confusion(:,:,[nTime:nTime+60]),3),4));
    figure
    plotConfMat(confusion_averagetime,{'red','green','blue'});
    string = sprintf('Alpha decode color: Timepoint between  %d and %d ms',timepoint_start,timepoint_end);
    axis square
    title(string,'FontSize', 14)
end