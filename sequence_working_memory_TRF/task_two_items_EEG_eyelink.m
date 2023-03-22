        sca;
        close all;
        clearvars;

        commandwindow;
        Screen('Preference', 'SkipSyncTests', 0);

        subnum = input('Input Sub No.: ');
        saveroot = ['data_sub_' num2str(subnum,'%02d') '/'];
        mkdir(saveroot);

        %===================open our screen====================
        sM = screenManager();
        sM.screen = max(Screen('Screens'));
        sM.pixelsPerCm = 27.5;
        sM.distance = 85.5;
        sM.debug = false;
        sM.blend = true;
        if sM.debug
            sM.bitDepth = '8bit';
        else
            sM.bitDepth = 'EnableBits++Bits++Output'; %EnableBits++Bits++Output EnableBits++Color++Output FloatingPoint32Bit
        end
        sM.backgroundColour = [0 0 0];
        sM.open; % OPEN THE SCREEN
        window = sM.win;

        waitframes = 1;
        ifi = Screen('GetFlipInterval', window);

        fprintf('\n--->>> Opened Screen %i : %s\n', sM.win, sM.fullName);

        %==============================setup eyelink==========================
        useEyeLink = 1;
        if useEyeLink == true
            eL = eyelinkManager('IP',[]);
            eL.name = ['data_sub_' num2str(subnum,'%02d')];
            fprintf('--->>> eL setup starting: %s\n', eL.name);
            eL.isDummy = false; %use dummy or real eyelink
            % eL.saveFile = [eL.name '.edf'];
            eL.recordData = true; %save EDF file
            eL.sampleRate = 500;
            eL.remoteCalibration = false; % manual calibration?
            eL.calibrationStyle = 'HV5'; % calibration style
            eL.modify.calibrationtargetcolour = [1 1 1];
            eL.modify.calibrationtargetsize = 0.5;
            eL.modify.calibrationtargetwidth = 0.05;
            eL.modify.waitformodereadytime = 500;
            eL.modify.devicenumber = -1; % -1 = use any keyboard
            %--------------fix parameters
            fixX = 0;
            fixY = 0;
            firstFixInit = 1;
            firstFixTime = 1;
            firstFixDiameter = 3;
            strictFixation = true;
            % X, Y, FixInitTime, FixTime, Radius, StrictFix
            updateFixationValues(eL, fixX, fixY, firstFixInit,...
                firstFixTime, firstFixDiameter, strictFixation);
            %sM.verbose = true; eL.verbose = true; sM.verbosityLevel = 10; eL.verbosityLevel = 4; %force lots of log output
            initialise(eL, sM); %use sM to pass screen values to eyelink
            setup(eL); % do setup and calibration
            fprintf('--->>> eL setup complete: %s\n', eL.name);
            WaitSecs('YieldSecs',0.25);
            getSample(eL); %make sure everything is in memory etc.
        else
            eL = [];
        end
        %-----initialise eyelink and draw fix spot




        %=======================set up the triggers===========================
        dPP = plusplusManager('sM',sM);
        dPP.strobeMode = 'eeg';
        dPP.open();
        dPP.verbose = true;

        %----------------------------------------------------------------------
        %                       Response Record
        %----------------------------------------------------------------------
        version = mod(subnum, 2);
        numBlock = 8;
        numTrial = 20;
        numTrialTotal = numBlock * numTrial;
        angledif = 30;   %Set the difference between angles of 2 cued items

        %----------------------------------------------------------------------
        %                       Blocklist Generator
        %----------------------------------------------------------------------
        fixchangelist_repeat = gen80; % 80 trials * 2 fix locations, 1-fix change
        fixchangelist_mirror = gen80;

        fixchangelist = {}; % 1*8 cells, one cell per block; 20 trials * 3 parameters (c1: fix1, c2: fix2, c3: 1-repeat 2-mirror) per cell
        idx_r = ones(10,1); idx_m = 2*idx_r; idx_rule = cat(1,idx_r,idx_m);
        for i = 1:numBlock
            blocklist = [cat(1,fixchangelist_repeat(10*(i-1)+1:i*10,:),fixchangelist_mirror(10*(i-1)+1:i*10,:)),idx_rule];
            while 1
                temp = randperm(20);
                idx_rule_shuffle = idx_rule(temp);
                idx_diff = [];
                for j = 1:18
                    idx_diff(j,1) = all(diff(idx_rule_shuffle(j:j+2)) == 0);
                end
                if all(idx_diff == 0)
                    blocklist = blocklist(temp,:);
                    fixchangelist{i} = blocklist;
                    break;
                end
            end
        end

        %----------------------------------------------------------------------
        %                       Color
        %----------------------------------------------------------------------
        dimColor = [0.6 0.6 0.6];
        feedbackColor = [0.5 0.5 0.5];
        cmplummax = hsv; % use hsv as default colormap
        cmplummax = rgb2hsv(cmplummax(1:63,:)); % turn rgb value into hsv

        %----------------------------------------------------
        %                       Color
        %----------------------------------------------------
        % root = randperm(21); root = root(1);
        % hsv h range: 0-180; red:0-10; yellow: 26-34; green: 35-77; blue: 100-124
        cmplumfirst = cmplummax;
        cmplumfirst(:,3) = cmplumfirst(:,3)*0.7;
        cmplumfirst = hsv2rgb(cmplumfirst); % turn hsv value into rgb
        root=1;
        c1 = cmplumfirst(root,:);
        c2 = cmplumfirst(root+21,:);
        c3 = cmplumfirst(root+42,:);
        colorlist = [c1; c2; c3];

        %----------------------------------------------------------------------
        %                      Keyboard Information
        %----------------------------------------------------------------------
        escape = KbName('ESCAPE');
        space = KbName('Space');

        %----------------------------------------------------------------------
        %                       Fixation
        %----------------------------------------------------------------------
        fixColor = sM.screenVals.white;
        baseFix = [0 0 10 10];
        maxFix = max(baseFix) * 1.01;
        centeredFix = CenterRectOnPointd(baseFix, sM.xCenter, sM.yCenter);

        %----------------------------------------------------------------------
        %                       3 Items & 2 Item Cues
        %----------------------------------------------------------------------
        %location of items
        numSides = 3;
        radius = 120;

        %----------------------------------------------------------------------
        %                       Probe
        %----------------------------------------------------------------------
        %draw probes
        baseProbe = [0 0 130 130];
        maxDiameter = max(baseProbe) * 1.01;


        %----------------------------------------------------------------------
        %                       Experimental Loop
        %----------------------------------------------------------------------
        % Beginning Instruction
        HideCursor(window);
        [~, secs, keyCode] = KbCheck;
        if version == 1
            imgStart = imread('Instruction_image/ImageInstructStar1.jpg');
        elseif version == 0
            imgStart = imread('Instruction_image/ImageInstructStar2.jpg');
        end
        imgStarTexture = Screen('MakeTexture', window, imgStart);
        Screen('DrawTexture', window, imgStarTexture, [], [], 0);
        Screen('Flip', window);
        KbStrokeWait;

        for block = 1:numBlock
            %----------------------------------------------------
            %                       Response Record
            %----------------------------------------------------
            blocklist = fixchangelist{1, block};
            item1locList = ones(numTrial, 4) *9;    %locations of 3 items
            item2locList = ones(numTrial, 4) *9;
            item3locList = ones(numTrial, 4) *9;
            angleList = ones(numTrial, 3) *9;
            probe1loc1List = ones(numTrial, 4) *9;   %locations of probes
            probe1loc2List = ones(numTrial, 4) *9;
            probe1loc3List = ones(numTrial, 4) *9;
            probe2loc1List = ones(numTrial, 4) *9;
            probe2loc2List = ones(numTrial, 4) *9;
            probe2loc3List = ones(numTrial, 4) *9;
            OT1locList = ones(numTrial, 4) *9;   %locations of 3 items in Order Test
            OT2locList = ones(numTrial, 4) *9;
            OT3locList = ones(numTrial, 4) *9;
            positionRecord = zeros(numTrial, 3);    %Order Test Response
            orderRecord = zeros(numTrial, 2) *9;
            timeRecord = ones(numTrial, 2) *9;
            timeRecord1 = ones(numTrial, 3) *9;
            correctResp = ones(numTrial, 2) *9;
            evaList = zeros(numTrial, 1);
            evaListall = zeros(numTrial, 2);
            fixEvaList = ones(numTrial, 2) *9;   %fixation change response
            fixResList = ones(numTrial, 2) *9;
            fixTimeRecord = ones(numTrial, 2) *9;

            Screen('TextSize', window, 40);
            Screen('TextFont', window, 'Courier');
            DrawFormattedText(window, ['Block_',num2str(block,'%02d'),'_start'] , 'center', 'center', sM.screenVals.white);
            Screen('Flip', window);
            WaitSecs(2);

            for trial = 1:numTrial
                % color luminance sequence
                for k = 1:6 % 3 item * 2 probe
                    lumr=(rand([1,600])*1); % 600 luminance sequence
                    alum=lumr-mean(lumr); % normalization
                    N=length(alum);
                    y=fft(alum,N);
                    y(find(real(y)==0))=1+0i;
                    % nl=0:N-1; % align to matlab index
                    % f=nl*100/N;% sequence of frequency
                    mag=abs(y);
                    realp=real(y);
                    imagp=imag(y);
                    weig=sqrt((repmat((mean(mag))^2,1,N))./(realp.^2+imagp.^2));% compute weight to normalize to mean. the case of mag=0 will lead error
                    unistim=(weig.*realp)+(weig.*imagp)*1i;
                    ulummid=ifft(unistim);
                    ulummid= ulummid+mean(lumr);
                    ulummid = (ulummid-min(ulummid))/(max(ulummid)- min(ulummid));
                    lumsequence(k,:) = ulummid;
                end
                cmplumlistdef = zeros([size(cmplummax),length(lumr)]); % create a list of 300 luminance-change value
                for k = 1:6
                    for cmpi = 1:length(lumr)
                        cmptmp = cmplummax;
                        cmptmp(:,3) = cmptmp(:,3)*lumsequence(k,cmpi);
                        cmptmp = hsv2rgb(cmptmp); % turn hsv value into rgb
                        cmplumlistdef(:,:,cmpi,k) = cmptmp;
                    end
                end
                %----------------------------------------------------
                %                       3 Items
                %----------------------------------------------------
                %color of items
                colornum = randperm(3);
                c1 = colorlist(colornum(1),:);
                c2 = colorlist(colornum(2),:);
                c3 = colorlist(colornum(3),:);
                c1tmp = rgb2hsv(c1);
                c1dim = [c1tmp(1), c1tmp(2), 1];
                c1dim = hsv2rgb(c1dim);
                c2tmp = rgb2hsv(c2);
                c2dim = [c2tmp(1), c2tmp(2), 1];
                c2dim = hsv2rgb(c2dim);
                c3tmp = rgb2hsv(c3);
                c3dim = [c3tmp(1), c3tmp(2), 1];
                c3dim = hsv2rgb(c3dim);
                dimcuelist = [c1dim; c2dim; c3dim];

                %location of items
                locDegree1 = randi(361) - 361;
                locDegree2 = locDegree1 + 360;
                anglesDegItemloc = linspace(locDegree1, locDegree2, numSides+1);
                anglesRadItemloc = anglesDegItemloc * (pi / 180);
                yPosVectorItemloc = sin(anglesRadItemloc) .* radius + sM.yCenter;
                xPosVectorItemloc = cos(anglesRadItemloc) .* radius + sM.xCenter;

                %draw items
                pic = ones(25,60,3)*255;
                angle1 = randi(180);
                angle2 = angle1 + angledif;
                angle3 = angle1 - angledif;
                angleList(trial,:) = [angle1, angle2, angle3];
                picn1 = pic_modify(pic,c1,angle1);
                picn2 = pic_modify(pic,c2,angle2);
                picn3 = pic_modify(pic,c3,angle3);
                item1loc = [xPosVectorItemloc(1)-size(picn1,2)/2 yPosVectorItemloc(1)-size(picn1,1)/2 xPosVectorItemloc(1)+size(picn1,2)/2 yPosVectorItemloc(1)+size(picn1,1)/2];
                item2loc = [xPosVectorItemloc(2)-size(picn2,2)/2 yPosVectorItemloc(2)-size(picn2,1)/2 xPosVectorItemloc(2)+size(picn2,2)/2 yPosVectorItemloc(2)+size(picn2,1)/2];
                item3loc = [xPosVectorItemloc(3)-size(picn3,2)/2 yPosVectorItemloc(3)-size(picn3,1)/2 xPosVectorItemloc(3)+size(picn3,2)/2 yPosVectorItemloc(3)+size(picn3,1)/2];
                item1locList(trial,:) = item1loc;
                item2locList(trial,:) = item2loc;
                item3locList(trial,:) = item3loc;
                trialitemloc = [item1loc; item2loc; item3loc];

                %location of item cues
                cueLoc = randperm(2);
                cueLoc1 = cueLoc(1);
                cueLoc2 = cueLoc(2);

                %cue pictures 1&2
                cue1piccolor = dimcuelist(cueLoc1,:);
                cue2piccolor = dimcuelist(cueLoc2,:);
                cue1picangle = angleList(trial,cueLoc1);
                cue2picangle = angleList(trial,cueLoc2);
                cue1picn = pic_modify(pic,cue1piccolor,cue1picangle);
                cue2picn = pic_modify(pic,cue2piccolor,cue2picangle);
                cue1picloc = trialitemloc(cueLoc1,:);
                cue2picloc = trialitemloc(cueLoc2,:);

                %draw item cues
                cueColor = sM.screenVals.white;
                radiusCue = 70;
                baseCue = [0 0 radiusCue*2 radiusCue*2];
                locCue1 = CenterRectOnPointd(baseCue, xPosVectorItemloc(cueLoc1), yPosVectorItemloc(cueLoc1));
                centerlocCue1 = [locCue1(1) + radiusCue, locCue1(2) + radiusCue];
                locCue2 = CenterRectOnPointd(baseCue, xPosVectorItemloc(cueLoc2), yPosVectorItemloc(cueLoc2));
                centerlocCue2 = [locCue2(1) + radiusCue, locCue2(2) + radiusCue];

                %-------------------------------------------------
                %                       Probe
                %-------------------------------------------------
                rotatenum = randperm(3);
                rotatedegree1 = 90 * rotatenum(1);
                rotatedegree2 = 90 * rotatenum(2);
                rotatedegree3 = 90 * rotatenum(3);
                %probe 1 locations
                anglesDegProbeloc1 = linspace(locDegree1-rotatedegree1, locDegree2-rotatedegree1, numSides+1);
                anglesRadProbeloc1 = anglesDegProbeloc1 * (pi / 180);
                yPosVectorProloc1 = sin(anglesRadProbeloc1) .* radius + sM.yCenter;
                xPosVectorProloc1 = cos(anglesRadProbeloc1) .* radius + sM.xCenter;

                %probe 2 locations
                anglesDegProbeloc2 = linspace(locDegree1-rotatedegree2, locDegree2-rotatedegree2, numSides+1);
                anglesRadProbeloc2 = anglesDegProbeloc2 * (pi / 180);
                yPosVectorProloc2 = sin(anglesRadProbeloc2) .* radius + sM.yCenter;
                xPosVectorProloc2 = cos(anglesRadProbeloc2) .* radius + sM.xCenter;

                %draw probe 1
                loc1Probe1 = CenterRectOnPointd(baseProbe, xPosVectorProloc1(1), yPosVectorProloc1(1));
                loc1Probe2 = CenterRectOnPointd(baseProbe, xPosVectorProloc1(2), yPosVectorProloc1(2));
                loc1Probe3 = CenterRectOnPointd(baseProbe, xPosVectorProloc1(3), yPosVectorProloc1(3));
                probe1loc1List(trial,:) = loc1Probe1;
                probe1loc2List(trial,:) = loc1Probe2;
                probe1loc3List(trial,:) = loc1Probe3;

                %draw probe 2
                loc2Probe1 = CenterRectOnPointd(baseProbe, xPosVectorProloc2(1), yPosVectorProloc2(1));
                loc2Probe2 = CenterRectOnPointd(baseProbe, xPosVectorProloc2(2), yPosVectorProloc2(2));
                loc2Probe3 = CenterRectOnPointd(baseProbe, xPosVectorProloc2(3), yPosVectorProloc2(3));
                probe2loc1List(trial,:) = loc2Probe1;
                probe2loc2List(trial,:) = loc2Probe2;
                probe2loc3List(trial,:) = loc2Probe3;

                %------Send strobe----%
                dPP.sendStrobe(248);  % trigger 248 stands for experiment begin
                Screen('Flip', window);

                %------Start------%
                Screen('FillOval', window, fixColor, centeredFix, maxFix);   %Fix
                Screen('PutImage', window, picn1, item1loc);   %3 Items
                Screen('PutImage', window, picn2, item2loc);
                Screen('PutImage', window, picn3, item3loc);
                Screen('Flip', window);
                WaitSecs(0.5);
                %------Item 1------%
                Screen('FillOval', window, fixColor, centeredFix, maxFix);   %Fix
                Screen('PutImage', window, picn1, item1loc);   %3 Items
                Screen('PutImage', window, picn2, item2loc);
                Screen('PutImage', window, picn3, item3loc);
                Screen('PutImage', window, cue1picn, cue1picloc);
                dashcircle(centerlocCue1,window,radiusCue,5,10);
                Screen('Flip', window);
                WaitSecs(1);
                %------Item 2------%
                Screen('FillOval', window, fixColor, centeredFix, maxFix); %Fix
                Screen('PutImage', window, picn1, item1loc);   %3 Items
                Screen('PutImage', window, picn2, item2loc);
                Screen('PutImage', window, picn3, item3loc);
                Screen('PutImage', window, cue2picn, cue2picloc);
                dashcircle(centerlocCue2,window,radiusCue,5,10);
                Screen('Flip', window);
                WaitSecs(1);
                %------Fixation 1------%

                if useEyeLink
                    resetFixation(eL);
                    trackerClearScreen(eL);
                    trackerDrawStimuli(eL);
                    trackerDrawFixation(eL); %draw fixation window on eyelink computer
                    edfMessage(eL,'V_RT MESSAGE END_FIX END_RT'); ... %this 3 lines set the trial info for the eyelink
                        edfMessage(eL,['TRIALID ' num2str(trial+20*(block-1)), num2str(1)]); ... %obj.getTaskIndex gives us which trial we're at
                        startRecording(eL);
                    statusMessage(eL,'INITIATE FIXATION...');
                    fixated = '';

                    Screen('FillOval', window, fixColor, centeredFix, maxFix); %Fix
                    Screen('DrawingFinished', sM.win); %tell PTB/GPU to draw
                    tFix = Screen('Flip',sM.win); %flip the buffer
                    WaitSecs(0.5)
                else
                    Screen('FillOval', window, fixColor, centeredFix, maxFix); %Fix
                    tFix = Screen('Flip',sM.win); %flip the buffer
                    WaitSecs(0.5);
                    fixated = 'fix';
                end

                %------Probe 1------%
                probeColorPos = randi(3);
                fixDimLastingTime = 0.2;
                fixDimStart = rand(1)*2.6 + 1;
                probe1cmplumlist1 = cmplumlistdef(:,:,:,1);
                probe1cmplumlist2 = cmplumlistdef(:,:,:,2);
                probe1cmplumlist3 = cmplumlistdef(:,:,:,3);

                %%%%% evaluate the luminance in sequence of color%%%%
                probe1lumc1 = lumsequence(1,:);
                probe1lumc2 = lumsequence(2,:);
                probe1lumc3 = lumsequence(3,:);
                probe1lum = [probe1lumc1; probe1lumc2; probe1lumc3];
                probe1lumorder01 = probe1lum(colornum(1),:);
                probe1lumorder02 = probe1lum(colornum(2),:);
                probe1lumorder03 = probe1lum(colornum(3),:);
                probe1lumorder = [probe1lumorder01; probe1lumorder02; probe1lumorder03];
                probe1lumsequence01 = probe1lumorder(cueLoc(1),:);
                probe1lumsequence02 = probe1lumorder(cueLoc(2),:);
                probe1lumsequence03 = probe1lumorder(3,:);
                probe1lumsequence = [probe1lumsequence01; probe1lumsequence02; probe1lumsequence03];
                probe1lumsequenceList(:,:,trial) = probe1lumsequence;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%

                %------Send strobe----%
                syncTime(eL);
                dPP.sendStrobe(249);  %trigger 249 stands for probe 1
                vbl = Screen('Flip', window);

                if probeColorPos == 1
                    timenow = GetSecs;
                    vbl(1) = timenow;
                    i = 1;
                    while vbl(i) < timenow + 5 && i<=600
                        probe1c1list = squeeze(probe1cmplumlist1(root,:,i));
                        probe1c2list = squeeze(probe1cmplumlist2(root+21,:,i));
                        probe1c3list = squeeze(probe1cmplumlist3(root+42,:,i));
                        Screen('FillOval', window, probe1c1list, loc1Probe1, maxDiameter);
                        Screen('FillOval', window, probe1c2list, loc1Probe2, maxDiameter);
                        Screen('FillOval', window, probe1c3list, loc1Probe3, maxDiameter);

                        %Fix
                        Screen('FillOval', window, fixColor, centeredFix, maxFix);
                        if blocklist(trial, 1) == 1
                            if GetSecs > timenow + fixDimStart && GetSecs < timenow + fixDimStart + fixDimLastingTime
                                Screen('FillOval', window, dimColor, centeredFix, maxFix);
                            end
                        end
                        i= i+1;
                        vbl(i) = Screen('Flip', window, vbl(i-1) + (waitframes - 0.5) * ifi);
                    end

                elseif probeColorPos == 2
                    timenow = GetSecs;
                    vbl(1) = timenow;
                    i=1;

                    while vbl(i) < timenow + 5 && i<=600
                        probe1c1list = squeeze(probe1cmplumlist1(root,:,i));
                        probe1c2list = squeeze(probe1cmplumlist2(root+21,:,i));
                        probe1c3list = squeeze(probe1cmplumlist3(root+42,:,i));
                        Screen('FillOval', window, probe1c2list, loc1Probe1, maxDiameter);
                        Screen('FillOval', window, probe1c3list, loc1Probe2, maxDiameter);
                        Screen('FillOval', window, probe1c1list, loc1Probe3, maxDiameter);

                        %Fix
                        Screen('FillOval', window, fixColor, centeredFix, maxFix);
                        if blocklist(trial, 1) == 1
                            if GetSecs > timenow + fixDimStart && GetSecs < timenow + fixDimStart + fixDimLastingTime
                                Screen('FillOval', window, dimColor, centeredFix, maxFix);
                            end
                        end
                        i = i+1;
                        vbl(i) = Screen('Flip', window, vbl(i-1) + (waitframes - 0.5) * ifi);
                    end
                else
                    timenow = GetSecs;
                    vbl(1) = timenow;
                    i=1;
                    while vbl(i) < timenow + 5 && i<=600
                        probe1c1list = squeeze(probe1cmplumlist1(root,:,i));
                        probe1c2list = squeeze(probe1cmplumlist2(root+21,:,i));
                        probe1c3list = squeeze(probe1cmplumlist3(root+42,:,i));
                        Screen('FillOval', window, probe1c3list, loc1Probe1, maxDiameter);
                        Screen('FillOval', window, probe1c1list, loc1Probe2, maxDiameter);
                        Screen('FillOval', window, probe1c2list, loc1Probe3, maxDiameter);

                        %Fix
                        Screen('FillOval', window, fixColor, centeredFix, maxFix);
                        if blocklist(trial, 1) == 1
                            if GetSecs > timenow + fixDimStart && GetSecs < timenow + fixDimStart + fixDimLastingTime
                                Screen('FillOval', window, dimColor, centeredFix, maxFix);
                            end
                        end
                        i = i+1;
                        vbl(i) = Screen('Flip', window, vbl(i-1) + (waitframes - 0.5) * ifi);
                    end
                end

                %------STOP RECORDING OF EYELINK
                if useEyeLink
                    statusMessage(eL,'End Probe 1');
                    edfMessage(eL,'MSG:endprobe1');
                    resetFixation(eL);
                    stopRecording(eL);
                    edfMessage(eL,'TRIAL_RESULT 1');
                    setOffline(eL);
                end
                %------------------------------

               
                %------Fix Test 1------%
                Screen('TextSize', window, 45);
                Screen('TextFont', window, 'Courier');
                DrawFormattedText(window, '?', 'center', 'center', sM.screenVals.white)
                Screen('Flip', window);
                fixEvaList(trial, 1) = 0;
                waitResp = 5;
                timenow = GetSecs;
                while GetSecs < timenow + waitResp
                    [~, secs, keyCode] = KbCheck;
                    [x, y, buttons] = GetMouse;
                    if version == 1
                        if buttons(1)
                            fprintf('left click 1\n');
                            fixResList(trial, 1) = 1;
                            fixTimeRecord(trial, 1) = GetSecs - timenow;
                            break
                        elseif buttons(3)
                            fprintf('right click 1\n');
                            fixResList(trial, 1) = 0;
                            fixTimeRecord(trial, 1) = GetSecs - timenow;
                            break
                        elseif keyCode(escape)
                            sca;
                            disp('*** Experiment terminated ***');
                            break
                        end
                    elseif version == 0
                        if buttons(3)
                            fprintf('right click 1\n');
                            fixResList(trial, 1) = 1;
                            fixTimeRecord(trial, 1) = GetSecs - timenow;
                            break
                        elseif buttons(1)
                            fprintf('left click 1\n');
                            fixResList(trial, 1) = 0;
                            fixTimeRecord(trial, 1) = GetSecs - timenow;
                            break
                        elseif keyCode(escape)
                            sca;
                            disp('*** Experiment terminated ***');
                            break
                        end
                    end
                end
                if fixResList(trial, 1) == blocklist(trial, 1)
                    fixEvaList(trial, 1) = 1;
                end

                %------Rule Cue------%
                if version == 1 || version ==0
                    if blocklist (trial, 3) == 1
                        Screen('TextSize', window, 40);
                        Screen('TextFont', window, 'Courier');
                        DrawFormattedText(window, 'repeat', 'center', 'center', sM.screenVals.white);
                        correctResp(trial, 1) = cueLoc1;
                        correctResp(trial, 2) = cueLoc2;
                    else
                        Screen('TextSize', window, 40);
                        Screen('TextFont', window, 'Courier');
                        DrawFormattedText(window, 'mirror', 'center', 'center', sM.screenVals.white);
                        correctResp(trial, 1) = cueLoc2;
                        correctResp(trial, 2) = cueLoc1;
                    end
                end
                Screen('Flip', window);
                WaitSecs(1);


                %------Fixation 2------%
                if useEyeLink
                    resetFixation(eL);
                    trackerClearScreen(eL);
                    trackerDrawStimuli(eL);
                    trackerDrawFixation(eL); %draw fixation window on eyelink computer
                    edfMessage(eL,'V_RT MESSAGE END_FIX END_RT'); ... %this 3 lines set the trial info for the eyelink
                        edfMessage(eL,['TRIALID ' num2str(trial+20*(block-1)), num2str(2)]); ... %obj.getTaskIndex gives us which trial we're at
                        startRecording(eL);
                    statusMessage(eL,'INITIATE FIXATION...');
                    fixated = '';

                    Screen('FillOval', window, fixColor, centeredFix, maxFix); %Fix
                    Screen('DrawingFinished', sM.win); %tell PTB/GPU to draw
                    tFix = Screen('Flip',sM.win); %flip the buffer
                    WaitSecs(0.5);
                else
                    Screen('FillOval', window, fixColor, centeredFix, maxFix); %Fix
                    tFix = Screen('Flip',sM.win); %flip the buffer
                    WaitSecs(0.5);
                    fixated = 'fix';
                end

                %------Probe 2------%
                probeColorPos = randi(3);
                fixDimLastingTime = 0.2;
                fixDimStart = rand(1)*2.6 + 1;
                probe2cmplumlist1 = cmplumlistdef(:,:,:,4);
                probe2cmplumlist2 = cmplumlistdef(:,:,:,5);
                probe2cmplumlist3 = cmplumlistdef(:,:,:,6);
                %%%%% evaluate the luminance in sequence of color%%%%
                probe2lumc1 = lumsequence(4,:);
                probe2lumc2 = lumsequence(5,:);
                probe2lumc3 = lumsequence(6,:);
                probe2lum = [probe2lumc1; probe2lumc2; probe2lumc3];
                probe2lumorder01 = probe2lum(colornum(1),:);
                probe2lumorder02 = probe2lum(colornum(2),:);
                probe2lumorder03 = probe2lum(colornum(3),:);
                probe2lumorder = [probe2lumorder01; probe2lumorder02; probe2lumorder03];
                probe2lumsequence01 = probe2lumorder(cueLoc(1),:);
                probe2lumsequence02 = probe2lumorder(cueLoc(2),:);
                probe2lumsequence03 = probe2lumorder(3,:);
                probe2lumsequence = [probe2lumsequence01; probe2lumsequence02; probe2lumsequence03];
                probe2lumsequenceList(:,:,trial) = probe2lumsequence;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%
                %------Send strobe----%
                syncTime(eL);
                dPP.sendStrobe(250); % trigger 250 stands for probe 2
                vbl = Screen('Flip', window);


                if probeColorPos == 1
                    timenow = GetSecs;
                    vbl(1) = timenow;
                    i=1;
                    while vbl(i) < timenow + 5 && i<=600
                        probe2c1list = squeeze(probe2cmplumlist1(root,:,i));
                        probe2c2list = squeeze(probe2cmplumlist2(root+21,:,i));
                        probe2c3list = squeeze(probe2cmplumlist3(root+42,:,i));
                        Screen('FillOval', window, probe2c1list, loc2Probe1, maxDiameter);
                        Screen('FillOval', window, probe2c2list, loc2Probe2, maxDiameter);
                        Screen('FillOval', window, probe2c3list, loc2Probe3, maxDiameter);

                        %Fix
                        Screen('FillOval', window, fixColor, centeredFix, maxFix);
                        if blocklist(trial, 2) == 1
                            if GetSecs > timenow + fixDimStart && GetSecs < timenow + fixDimStart + fixDimLastingTime
                                Screen('FillOval', window, dimColor, centeredFix, maxFix);
                            end
                        end
                        i = i + 1;
                        vbl(i) = Screen('Flip', window, vbl(i-1) + (waitframes - 0.5) * ifi);
                    end

                elseif probeColorPos == 2
                    timenow = GetSecs;
                    vbl(1) = timenow;
                    i=1;
                    while vbl(i) < timenow + 5 && i<=600
                        probe2c1list = squeeze(probe2cmplumlist1(root,:,i));
                        probe2c2list = squeeze(probe2cmplumlist2(root+21,:,i));
                        probe2c3list = squeeze(probe2cmplumlist3(root+42,:,i));
                        Screen('FillOval', window, probe2c2list, loc2Probe1, maxDiameter);
                        Screen('FillOval', window, probe2c3list, loc2Probe2, maxDiameter);
                        Screen('FillOval', window, probe2c1list, loc2Probe3, maxDiameter);

                        %Fix
                        Screen('FillOval', window, fixColor, centeredFix, maxFix);
                        if blocklist(trial, 2) == 1
                            if GetSecs > timenow + fixDimStart && GetSecs < timenow + fixDimStart + fixDimLastingTime
                                Screen('FillOval', window, dimColor, centeredFix, maxFix);
                            end
                        end
                        i = i+1;
                        vbl(i) = Screen('Flip', window, vbl(i-1) + (waitframes - 0.5) * ifi);
                    end
                else
                    timenow = GetSecs;
                    vbl(1) = timenow;
                    i=1;
                    while vbl(i) < timenow + 5 && i<=600
                        probe2c1list = squeeze(probe2cmplumlist1(root,:,i));
                        probe2c2list = squeeze(probe2cmplumlist2(root+21,:,i));
                        probe2c3list = squeeze(probe2cmplumlist3(root+42,:,i));
                        Screen('FillOval', window, probe2c3list, loc2Probe1, maxDiameter);
                        Screen('FillOval', window, probe2c1list, loc2Probe2, maxDiameter);
                        Screen('FillOval', window, probe2c2list, loc2Probe3, maxDiameter);

                        %Fix
                        Screen('FillOval', window, fixColor, centeredFix, maxFix);
                        if blocklist(trial, 2) == 1
                            if GetSecs > timenow + fixDimStart && GetSecs < timenow + fixDimStart + fixDimLastingTime
                                Screen('FillOval', window, dimColor, centeredFix, maxFix);
                            end
                        end
                        i = i+1;
                        vbl(i) = Screen('Flip', window, vbl(i-1) + (waitframes - 0.5) * ifi);
                    end
                end

                %------STOP RECORDING OF EYELINK
                if useEyeLink
                    statusMessage(eL,'End Probe 2');
                    edfMessage(eL,'MSG:endprobe2');
                    resetFixation(eL);
                    stopRecording(eL);
                    edfMessage(eL,'TRIAL_RESULT 1');
                    setOffline(eL);
                end
                %------------------------------

                

                %------Fix Test 2------%
                Screen('TextSize', window, 45);
                Screen('TextFont', window, 'Courier');
                DrawFormattedText(window, '?', 'center', 'center', sM.screenVals.white);
                Screen('Flip', window);
                fixEvaList(trial, 2) = 0;
                waitResp = 5;
                timenow = GetSecs;
                while GetSecs < timenow + waitResp
                    [~, secs, keyCode] = KbCheck;
                    [x, y, buttons] = GetMouse;
                    if version == 1
                        if buttons(1)
                            fprintf('left click 2\n');
                            fixResList(trial, 2) = 1;
                            fixTimeRecord(trial, 2) = GetSecs - timenow;
                            break
                        elseif buttons(3)
                            fprintf('right click 2\n');
                            fixResList(trial, 2) = 0;
                            fixTimeRecord(trial, 2) = GetSecs - timenow;
                            break
                        elseif keyCode(escape)
                            sca;
                            disp('*** Experiment terminated ***');
                            break
                        end
                    elseif version == 0
                        if buttons(3)
                            fprintf('right click 2\n');
                            fixResList(trial, 2) = 1;
                            fixTimeRecord(trial, 2) = GetSecs - timenow;
                            break
                        elseif buttons(1)
                            fprintf('left click 2\n');
                            fixResList(trial, 2) = 0;
                            fixTimeRecord(trial, 2) = GetSecs - timenow;
                            break
                        elseif keyCode(escape)
                            sca;
                            disp('*** Experiment terminated ***');
                            break
                        end
                    end
                end
                if fixResList(trial, 2) == blocklist(trial, 2)
                    fixEvaList(trial, 2) = 1;
                end

                %------Order Test------%
                anglesDegTestloc = linspace(locDegree1-rotatedegree3, locDegree2-rotatedegree3, numSides+1);
                anglesRadTestloc = anglesDegTestloc * (pi / 180);
                yPosVectorTestloc = sin(anglesRadTestloc) .* radius + sM.yCenter;
                xPosVectorTestloc = cos(anglesRadTestloc) .* radius + sM.xCenter;

                positionnum = randperm(3);
                item1loc = [xPosVectorTestloc(positionnum(1))-size(picn1,2)/2 yPosVectorTestloc(positionnum(1))-size(picn1,1)/2 xPosVectorTestloc(positionnum(1))+size(picn1,2)/2 yPosVectorTestloc(positionnum(1))+size(picn1,1)/2];
                item2loc = [xPosVectorTestloc(positionnum(2))-size(picn2,2)/2 yPosVectorTestloc(positionnum(2))-size(picn2,1)/2 xPosVectorTestloc(positionnum(2))+size(picn2,2)/2 yPosVectorTestloc(positionnum(2))+size(picn2,1)/2];
                item3loc = [xPosVectorTestloc(positionnum(3))-size(picn3,2)/2 yPosVectorTestloc(positionnum(3))-size(picn3,1)/2 xPosVectorTestloc(positionnum(3))+size(picn3,2)/2 yPosVectorTestloc(positionnum(3))+size(picn3,1)/2];
                OT1locList(trial,:) = item1loc;
                OT2locList(trial,:) = item2loc;
                OT3locList(trial,:) = item3loc;
                xPosVectorItem1 = [item1loc(1), item1loc(3)];
                yPosVectorItem1 = [item1loc(2), item1loc(4)];
                xPosVectorItem2 = [item2loc(1), item2loc(3)];
                yPosVectorItem2 = [item2loc(2), item2loc(4)];
                xPosVectorItem3 = [item3loc(1), item3loc(3)];
                yPosVectorItem3 = [item3loc(2), item3loc(4)];

                picn1 = pic_modify(pic,sM.screenVals.white,angle1);
                picn2 = pic_modify(pic,sM.screenVals.white,angle2);
                picn3 = pic_modify(pic,sM.screenVals.white,angle3);
                SetMouse(sM.xCenter, sM.yCenter, window);
                ShowCursor('Arrow', window);
                Screen('FillOval', window, fixColor, centeredFix, maxFix); %Fix
                Screen('PutImage', window, picn1, item1loc);   %3 Items
                Screen('PutImage', window, picn2, item2loc);
                Screen('PutImage', window, picn3, item3loc);
                Screen('Flip', window,[],1);
                order = 0;
                in1 = 0;
                in2 = 0;
                in3 = 0;
                count = 0;
                waitResp = 7;
                timenow = GetSecs;

                while GetSecs < timenow + waitResp
                    [x, y, buttons] = GetMouse;
                    in1 = inpolygon(x, y, xPosVectorItem1, yPosVectorItem1);
                    in2 = inpolygon(x, y, xPosVectorItem2, yPosVectorItem2);
                    in3 = inpolygon(x, y, xPosVectorItem3, yPosVectorItem3);

                    if sum(buttons) > 0
                        fprintf('1st click\n');
                        count = count + 1;
                        if in1 == 1 && count == 1
                            picn1 = pic_modify(pic,feedbackColor,angle1);
                            Screen('PutImage', window, picn1, item1loc);
                            Screen('Flip', window,[],1);
                            positionRecord(trial, 1) = 1;
                            order = sum(positionRecord(trial,:));
                            orderRecord(trial, order) = 1;
                            timeRecord(trial, order) = GetSecs - timenow;
                            timeRecord1(trial, 1) = GetSecs - timenow;
                            if order == 1
                                fprintf('hellA1\n');
                            else
                                fprintf('hellB1\n');
                            end
                        elseif in2 == 1 && count == 1
                            picn2 = pic_modify(pic,feedbackColor,angle2);
                            Screen('PutImage', window, picn2, item2loc);
                            Screen('Flip', window,[],1);
                            positionRecord(trial, 2) = 1;
                            order = sum(positionRecord(trial,:));
                            orderRecord(trial, order) = 2;
                            timeRecord(trial, order) = GetSecs - timenow;
                            timeRecord1(trial, 2) = GetSecs - timenow;
                            if order == 1
                                fprintf('hellA2\n');
                            else
                                fprintf('hellB2\n');
                            end
                        elseif in3 == 1 && count == 1
                            picn3 = pic_modify(pic,feedbackColor,angle3);
                            Screen('PutImage', window, picn3, item3loc);
                            Screen('Flip', window,[],1);
                            positionRecord(trial, 3) = 1;
                            order = sum(positionRecord(trial,:));
                            orderRecord(trial, order) = 3;
                            timeRecord(trial, order) = GetSecs - timenow;
                            timeRecord1(trial, 3) = GetSecs - timenow;
                            if order == 1
                                fprintf('hellA3\n');
                            else
                                fprintf('hellB3\n');
                            end
                        else
                            while sum(buttons) > 0
                                in1 = 0;
                                in2 = 0;
                                in3 = 0;
                                count = 0;
                                [x, y, buttons] = GetMouse;
                            end
                            continue
                        end

                        while sum(buttons) > 0 && order == 1
                            in1 = 0;
                            in2 = 0;
                            in3 = 0;
                            count = 0;
                            [x, y, buttons] = GetMouse;
                        end

                        if order == 2
                            break
                        end

                    end
                end
                % ------------------------------------

                fprintf('end check\n');
                %Evaluate the response
                if orderRecord(trial, 1) == correctResp(trial, 1) && orderRecord(trial, 2) == correctResp(trial, 2)
                    evaList(trial) = 1;
                    dPP.sendStrobe(251);  %trigger 251 stands for right order response
                    Screen('Flip', window);
                else
                    evaList(trial) = 0;
                    dPP.sendStrobe(252);  %trigger 252 stands for wrong order response
                    Screen('Flip', window);
                end


                %evaluate all the response

                for i = 1:2
                    if orderRecord(trial, i) == correctResp(trial, i)
                        evaListall(trial, i) = 1;
                    else
                        evaListall(trial, i) = 0;
                    end
                end

                %ITI
                HideCursor(window);
                Screen('Flip', window);
                Screen('FillOval', window, fixColor, centeredFix, maxFix); %Fix
                Screen('Flip', window);
                itiLength = 2 + rand(1);
                timenow = GetSecs;
                while GetSecs < timenow + itiLength
                    [~, secs, keyCode] = KbCheck;
                    if keyCode(escape)
                        sca;
                        disp('*** Experiment terminated ***');
                        break
                    end
                end
            end

            %Break between blocks
            subfile = ['Sub_' num2str(subnum,'%02d') '_' date '_block' num2str(block,'%02d')];
            save([saveroot subfile],'','');
			

			
            if block<numBlock
                waitResp = 120;
                timenow = GetSecs;
                imgBreak = imread('Instruction_image/ImageInstructBreak.jpg');
                imgBreakTexture = Screen('MakeTexture', window, imgBreak);
                Screen('DrawTexture', window, imgBreakTexture, [], [], 0);
                Screen('Flip', window);
                while GetSecs < timenow + waitResp
                    [~, secs, keyCode] = KbCheck;
                    if keyCode(space)
                        break
                    end
                end
            end
		end
		
		if useEyeLink == true
			subfileall = ['Sub_' num2str(subnum,'%02d') '_' date];
				eL.saveFile = [saveroot subfileall '.edf'];
				close(eL);
			end

        imgEnd = imread('Instruction_image/ImageInstructEnd.jpg');
        imgEndTexture = Screen('MakeTexture', window, imgEnd);
        Screen('DrawTexture', window, imgEndTexture, [], [], 0);
        Screen('Flip', window);
        KbStrokeWait;
        sca;


        %----------------------------------------------------------------------
        %                       Functions
        %----------------------------------------------------------------------
        function picn = pic_modify(pic,color,angle)

        picn = zeros(size(pic));
        for x = 1:size(pic,1)
            for y = 1:size(pic,2)
                c = pic(x,y,:);
                if any(c)
                    picn(x,y,:) = color;
                end
            end
        end

        picn = imrotate(picn,angle);
        end

        function fixchangelist = gen80

        while 1
            fixchangelist_ori = [ones(20,1);zeros(60,1)];
            while 1
                fixchangelist = [fixchangelist_ori(randperm(80)),fixchangelist_ori(randperm(80))];
                for i = 1:size(fixchangelist,1)
                    list(i,1) = all(fixchangelist(i,:));
                end
                if length(find(list==1)) == 10
                    break
                end
            end
            idx_both = find(list == 1);
            idx_1st = setdiff(find(fixchangelist(:,1) == 1),idx_both);
            idx_2nd = setdiff(find(fixchangelist(:,2) == 1),idx_both);
            fixchangelist = zeros(80,3);
            fixchangelist(idx_1st,1) = 1;
            fixchangelist(idx_2nd,2) = 1;
            fixchangelist(idx_both,3) = 1;
            for i = 1:8
                temp = fixchangelist((i-1)*10+1:i*10,:);
                summ(i) = sum(temp(:));
            end
            if all(summ <= 5) && all(summ >= 2)
                break;
            end
        end
        fixchangelist(fixchangelist(:,3) == 1,1:2) = 1;
        fixchangelist = fixchangelist(:,1:2);

        end

        function dashcircle(imgcenter,windowPtr,radius,width,angle)
        for i = 1 : round(360/angle)
            Screen('FrameArc',windowPtr,[1,1,1],[imgcenter(1)-radius,imgcenter(2)-radius,imgcenter(1)+radius,imgcenter(2)+radius],angle*(2*(i-1)),angle,width) ;
        end
        end
