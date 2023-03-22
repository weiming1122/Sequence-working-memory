clc;
clear;
close all;

subList = [1];
Nsub = length(subList);

Nblock = 3; % cross-validation
Nitr = 10; % iteration
Ntp = 441; % # of time points
NBins = 8; % # of location bin

tm = -200:5:2000;

%create empty matrix
AverageAccuracy = zeros(Nsub,Ntp);

for sub = 1:Nsub
        DecodingAccuracy = zeros(Ntp,Nblock,Nitr);
        %% load SVM_ECOC output files
        fileLocation = pwd;
        readThis =strcat(fileLocation,'/Orientation_Results_Alphabased_',num2str(subList(sub),'%02d'),'.mat');
        load(readThis);
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
     AverageAccuracy(sub,:) =smoothed; % average across iteration and block
     
end %End of subject

%compute average accuracy across participants and SE of the mean.
subAverage = squeeze(mean(AverageAccuracy,1)); 

%% plot
figure(1)
plot(tm,mean(AverageAccuracy,1),'LineWidth',2,'color','k');
xlim([-200 2000])
xlabel('Time (ms)','fontsize',14);ylabel('Decoding Accuracy','fontsize',14)
hold on
h = line(tm,1/8* ones(1,Ntp));
h.LineStyle = '--';
h.Color = [0.1,0.1,0.1];
set(gca,'linewidth',1,'fontsize',13);
title('Alpha-based decoding: Orientation','fontsize',14)
hold off
saveas(figure(1),'decoding_orientation_Alpha','png')