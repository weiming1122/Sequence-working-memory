    function sequence_memory_3items_part01
    % Sequence memory of orientations of 3 items
    close all;
    clc;
    clear all;

    warning('off','MATLAB:dispatcher:InexactMatch');                           % turn off case-mismatch manager (it's annoying)
    KbName('UnifyKeyNames');

    % Build a GUI for subject information
    prompt = {'Subject Number', 'Random Seed'};                                % What info do we want from the subject?
    s = ClockRandSeed;                                                         % grab a seed
    defAns = {'',num2str(s.Seed)};                                           % fill in some default answers
    % defAns = {'',num2str(s)};                                                  % fill in some default answers                                          % fill in some default answers
    box = inputdlg(prompt, 'Enter Subject Information', 1, defAns);            % actually make GUI
    if length(box) == length(defAns)                                           % simple check for enough input, otherwise bail
        p.subNum = num2str(box{1});  p.rndSeed = str2num(box{2});
        rand('state',p.rndSeed);                                               % seed the generator
    else
        return;
    end

    % Determine if Practice
    prompt = {'Practice (1=yes,0=no)'};
    defAns = {''};                                                              % fill in some default answers
    box = inputdlg(prompt, 'Enter Subject Information', 1, defAns);             % actually make GUI
    if length(box) == length(defAns)                                            % simple check for enough input, otherwise bail
        p.practice = str2num(box{1});                                           % seed the generator
    else
        return;
    end

    % Turn of keyboard echoing
    ListenChar(2);

    %--------------------------------------------------------------------------
    % General Experimental Parameters
    %--------------------------------------------------------------------------
    p.portCode = 1;           % (1: EEG recording; 0: no EEG recording)
    p.eyeLinkCode = 0;        % (1: eyeLink recording; 0: no eyeLink recording)

    p.nBlocks = 5;                                     % number of blocks
    p.nColors = 3;                                     % number of colors
    p.nOrientations = 12;                              % number of orientations
    p.nReps = 2;                                       % number of repetitions per color and per orientation for each block
    p.nTrials = p.nColors*p.nOrientations*p.nReps;     % number of trials per block

    % Dimensions for Stims
    p.vDist = 85;                                     % viewing distance (cm)
    p.px = 0.0248;                                     % pixel size (cm) for monitors: 1920*1080, 21.5 inch
    p.rimSize = deg2pix(6,p.vDist,p.px);               % diameter of rim
    p.rimThick = deg2pix(0.2,p.vDist,p.px);            % thickness of rim
    p.fixSize = deg2pix(0.2,p.vDist,p.px);             % fixation size
    p.triangleH1 = deg2pix(3,p.vDist,p.px);         % height of triangle part of stimulus sample
    p.triangleH2 = deg2pix(0.8,p.vDist,p.px);         % bottom of triangle part of stimulus sample
    p.semicircle = deg2pix(0.8,p.vDist,p.px);          % diameter of semicircle part of stimulus sample


    % Stimulus orientation for forward encoding model (FEM) analysis
    p.nChans = p.nOrientations;
    p.orienChannels = linspace(15,375-360/p.nChans,p.nChans);         % orientation vector: 15:30:345

    % Timing
    p.refreshRate = 60;
    p.refreshCycle = 1/p.refreshRate;
    p.ITI = [1.0:p.refreshCycle:1.5];                % ITI ranges between 1000 to 1500 ms
    p.fixDur = 0.5;                                  % duration of pre-sample fixation
    p.stimDur = 0.2;                                 % duration of sample stimuli
    p.delay = 1.3;                                   % duration of delay
    p.waitOrienResp = 4;                             % duration of waiting for response orientation

    % Color information
    p.black = [0 0 0];                     % color of black
    p.white = [215 215 215];               % color of white
    p.red = [215 0 0];                     % color of red
    p.blue = [0 0 215];                    % color of blue
    p.green = [0 215 0];                   % color of green

    p.foreCol = p.black;                   % foreground window color is black
    p.txtCol = p.white;                    % text color is white
    p.fixCol = p.white;                    % fixation color is white
    p.rimCol = p.white;                    % rim color is white

    if p.portCode == 1
        % trigger number
        p.itiTrigger = 240;
        p.fixTrigger = 241;
        p.stimTrigger = 242;
        p.delayTrigger = 243;
        p.ringTrigger = 244;
        p.orinRepTrigger = 245;
        p.triggerReset = 0;
    end

    %--------------------------------------------------------------------------
    % Build the stimulus display
    %--------------------------------------------------------------------------
    AssertOpenGL;                                                   % Make sure we're on an OpenGL-compatible machine
    s = max(Screen('Screens'));                                     % Grab a handle for the display ID (should return 0 for single-monitor setups)
    Screen('Preference', 'SkipSyncTests', 0);
    [w,p.sRect] = Screen('OpenWindow',s,p.foreCol);                 % full-screen
    p.ifi = Screen('GetFlipInterval', w);                           % estimate the monitor flip interval
    HideCursor;
    Priority(MaxPriority(w));                                       % set priority to max to discourage interruptions

    % Turn on alpha blending. this makes drawing the stims much easier.
    % Type Screen('BlendingFunction?') for more info.
    Screen('BlendFunction', w, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

    % Compute and store the center of the screen
    p.xCenter = (p.sRect(3)-p.sRect(1))/2;
    p.yCenter = (p.sRect(4)-p.sRect(2))/2;

    % Foreground rectangle
    p.foreRect = p.sRect;                    % set the foreground/foreRect to the full screen/p.sRect

    % Center fixation rectangle on the main display (foreRect)
    p.fixRect = CenterRect([0 0 p.fixSize p.fixSize],p.foreRect);

    % ring position
    p.ringPos = [p.xCenter - p.rimSize/2, p.yCenter - p.rimSize/2, ...
        p.xCenter + p.rimSize/2, p.yCenter + p.rimSize/2]';

    % make a dummy call to GetSecs to load the .dll before we need it
    dummy = GetSecs; clear dummy;

    %--------------------------------------------------------------------------
    % eeg trigger setup
    %--------------------------------------------------------------------------
    if p.portCode == 1
        ioObj = io64; % initialize inpoutx64.dll device driver
        if ( io64(ioObj) == 0 )
            disp(' ');
            disp('inpoutx64.dll successfully installed.')
        else
            disp('inpoutx64.dll installation failed.')
            return
        end
        address = hex2dec('D000');   % LPT3 Status port addr
        x=io64(ioObj,address);       % installation the io64 Driver
    end

    %--------------------------------------------------------------------------
    % setup eyelink
    %--------------------------------------------------------------------------
    if p.eyeLinkCode == 1
        % need temp name because Eyelink only know hows to save names with 8 chars or less.
        % Will change name using matlab's moveFile later.
        fileName = ['sub',num2str(p.subNum,'%02d'),'_part01'];
        savePath = fullfile(pwd,'eyelink_data');
        mkdir(savePath);
        tempName = 'TEMP1';
        dummymode=0;

        eL=EyelinkInitDefaults(w);

        % Initialization of the connection with the Eyelink Tracker.
        % exit program if this fails.
        if ~EyelinkInit(dummymode)
            fprintf('Eyelink Init aborted.\n');
            cleanup;  % cleanup function
            Eyelink('ShutDown');
            Screen('CloseAll');
            return;
        end

        i = Eyelink('Openfile', tempName);

        if i~=0
            fprintf('Cannot create EDF file ''%s'' ', fileName);
            cleanup;
            Eyelink('ShutDown');
            Screen('CloseAll');
            return
        end

        %   SET UP TRACKER CONFIGURATION
        Eyelink('command', 'calibration_type = HV9');
        Eyelink('command', 'online_dcorr_refposn = %1d, %1d', p.xCenter, p.yCenter);
        Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,AREA');


        % you must call this function to apply the changes from above
        EyelinkUpdateDefaults(eL);

        % Calibrate the eye tracker
        EyelinkDoTrackerSetup(eL);

        % do a final check of calibration using driftcorrection
        EyelinkDoDriftCorrection(eL);

        Eyelink('StartRecording');

        Eyelink('message', 'SYNCTIME');	 	     % zero-plot time for EDFVIEW
        eye_used = Eyelink('EyeAvailable');      % get eye that's tracked

        if eye_used == eL.BINOCULAR              % if both eyes are tracked
            eye_used = eL.LEFTEYE;               % use left eye
        end

        errorCheck=Eyelink('checkrecording'); 		% Check recording status

        if(errorCheck~=0)
            fprintf('Eyelink checked wrong status.\n');
            cleanup;  % cleanup function
            Eyelink('ShutDown');
            Screen('CloseAll');
        end
    end
    %--------------------------------------------------------------------------
    % Begin block loop
    %--------------------------------------------------------------------------

    for b = 1:p.nBlocks

        % Build an output directory & check to make sure it doesn't already exist

        p.root = pwd;

        if ~exist([p.root, '\Subject Data\'], 'dir')
            mkdir([p.root, '\Subject Data\']);
        end

        % Build an output file and check to make sure that it doesn't exist yet either
        fName = [p.root, '\Subject Data\', 'sub',num2str(p.subNum,'%02d'), '_OrienWM_part01_', num2str(b,'%02d'), '.mat'];

        if exist(fName)
            Screen('CloseAll');
            msgbox('File already exists!', 'modal');
            return;
        end

        %----------------------------------------------------------------------
        % Timing vectors
        %----------------------------------------------------------------------

        time.iti.vblstamp = nan(1,p.nTrials);
        time.iti.onset = nan(1,p.nTrials);
        time.iti.flipstamp = nan(1,p.nTrials);
        time.iti.missed = nan(1,p.nTrials);

        time.fixation.vblstamp = nan(1,p.nTrials);
        time.fixation.onset = nan(1,p.nTrials);
        time.fixation.flipstamp = nan(1,p.nTrials);
        time.fixation.missed = nan(1,p.nTrials);

        time.stim.vblstamp = nan(1,p.nTrials);
        time.stim.onset = nan(1,p.nTrials);
        time.stim.flipstamp = nan(1,p.nTrials);
        time.stim.missed = nan(1,p.nTrials);

        time.delay.vblstamp = nan(1,p.nTrials);
        time.delay.onset = nan(1,p.nTrials);
        time.delay.flipstamp = nan(1,p.nTrials);
        time.delay.missed = nan(1,p.nTrials);

        time.orientationResp.vblstamp = nan(1,p.nTrials);
        time.orientationResp.onset = nan(1,p.nTrials);
        time.orientationResp.flipstamp = nan(1,p.nTrials);
        time.orientationResp.missed = nan(1,p.nTrials);

        %----------------------------------------------------------------------
        % Control parameters
        %----------------------------------------------------------------------

        % preallocation stimuli label vectors
        stim.orienChanLab = nan(1,p.nTrials);               % stimulus orientation label
        stim.angle = nan(p.nTrials,1);                      % stimulus orientation angle
        stim.Degree = nan(p.nTrials,1);                     % stimulus orientation degree
        stim.A1Ordin =  nan(p.nTrials,2);                   % ordination of first vertex of triangle (part of stimulus)
        stim.A2Ordin =  nan(p.nTrials,2);                   % ordination of second vertex of triangle (part of stimulus)
        stim.A3Ordin =  nan(p.nTrials,2);                   % ordination of third vertex of triangle (part of stimulus)
        stim.triangleOrdin = nan(3,2,p.nTrials);            % ordination of triangle (part of stimulus)
        stim.arcRectOrdin = nan(p.nTrials,4);               % ordination of rect of semicircle (part of stimulus)
        stim.color = nan(1,p.nTrials);                      % stimulus color label
        stim.pColor = nan(p.nTrials,3);                     % stimulus color

        % preallocation orientation response vectors
        stim.tAngle = nan(1,p.nTrials);                     % target value (same as stim.probedPos)
        stim.rAngle = nan(1,p.nTrials);                     % actual reported value
        stim.rt = nan(1,p.nTrials);                         % response time
        stim.rLoc = nan(p.nTrials,2);                       % coordinates of mouseclick
        stim.angleOffset = nan(1,p.nTrials);                % offset error

        % scramble trial order
        stim.rndInd = randperm(p.nTrials);

        % create orientation Channel sequence and color sequence for stimuli
        cond.red = p.nChans;
        cond.green = p.nChans;
        cond.blue = p.nChans;

        [tStructure] = buildTrialStructure(cond,p);

        stim.tStructure = tStructure(:,stim.rndInd);
        stim.orienChanLab = tStructure(1,stim.rndInd);           % 1:1:12 (15:30:345)
        stim.color = tStructure(2,stim.rndInd);               % 1:red; 2:green; 3:blue

        clear cond

        %----------------------------------------------------------------------
        % Define Stimulus Parameters Before starting block
        %----------------------------------------------------------------------

        cnt=1;

        for t = 1:p.nTrials

            if rem(cnt,p.nTrials) == 0

                % Display for Block loading

                text1 = ['Block ' num2str(b) ' Loading'];
                tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, text1))/2 p.yCenter-120];
                text2 = [num2str(round2(cnt/p.nTrials,1e-2)*100), ' %'];
                tCenter2 = [p.xCenter-RectWidth(Screen('TextBounds', w, text2))/2 p.yCenter-70];

                Screen('FillRect',w,p.foreCol,p.foreRect);             % Draw the foreground window
                Screen('FillOval',w,p.fixCol,p.fixRect);               % Draw the fixation point
                Screen('DrawText', w, text1, tCenter1(1), tCenter1(2), p.txtCol);
                Screen('DrawText', w, text2, tCenter2(1), tCenter2(2), p.txtCol);
                Screen('DrawingFinished', w);
                Screen('Flip', w);

                [keyIsDown,secs,keyCode]=KbCheck;

                if keyIsDown
                    kp = find(keyCode); kp=kp(1);
                    if kp == 27                  % 'Esc'
                        Screen('CloseAll')
                        return
                    end
                end
            end

            % set ITI for trial
            stim.ITI(t) = randsample(p.ITI,1);

            % define location of stimuli
            stim.angle(t,:) = p.orienChannels(stim.orienChanLab(t));
            stim.Degree(t,:) = stim.angle(t,:)*pi/180;
            stim.A1Ordin(t,:) = [(p.xCenter + p.triangleH1 * cos(stim.Degree(t,:))),(p.yCenter + p.triangleH1 * sin(stim.Degree(t,:)))];
            stim.A2Ordin(t,:) = [(p.xCenter + p.triangleH2 * cos(stim.Degree(t,:) + 90*pi/180)),(p.yCenter + p.triangleH2 * sin(stim.Degree(t,:) + 90*pi/180))];
            stim.A3Ordin(t,:) = [(p.xCenter + p.triangleH2 * cos(stim.Degree(t,:) - 90*pi/180)),(p.yCenter + p.triangleH2 * sin(stim.Degree(t,:) - 90*pi/180))];
            stim.triangleOrdin(:,:,t) = [ stim.A1Ordin(t,:); stim.A2Ordin(t,:); stim.A3Ordin(t,:)];
            stim.arcRectOrdin(t,:) = [(p.xCenter - p.semicircle), (p.yCenter - p.semicircle), (p.xCenter + p.semicircle),(p.yCenter + p.semicircle)];

            % define colors of stimuli
            if stim.color(t) == 1
                stim.pColor(t,:) = p.red;    % 1:red
            elseif stim.color(t) == 2
                stim.pColor(t,:) = p.green;  % 2:green
            else
                stim.pColor(t,:) = p.blue;   % 3:blue
            end

            cnt=cnt+1;
        end

        % Display for Block loading
        text1 = ['You have ' num2str(p.nBlocks-(b-1)) ' blocks remaining.'];
        tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, text1))/2 p.yCenter-120];
        text2 = 'Press the spacebar to continue';
        tCenter2 = [p.xCenter-RectWidth(Screen('TextBounds', w, text2))/2 p.yCenter-70];
        Screen('FillRect',w,p.foreCol,p.foreRect);             % Draw the foreground window
        Screen('FillOval',w,p.fixCol,p.fixRect);               % Draw the fixation point
        Screen('DrawText', w, text1, tCenter1(1), tCenter1(2), p.txtCol);
        Screen('DrawText', w, text2, tCenter2(1), tCenter2(2), p.txtCol);
        Screen('DrawingFinished', w);
        Screen('Flip', w);
        WaitSecs(.5);

        % Wait for a spacebar press to continue
        while 1
            [keyIsDown,secs,keyCode]=KbCheck;
            if keyIsDown
                kp = find(keyCode);
                if kp == 32              % 'Spacebar'
                    break;
                end
            end
        end

        %----------------------------------------------------------------------
        % Begin Trial Loop
        %----------------------------------------------------------------------

        for t = 1:p.nTrials

            %----------------------------------------------
            % ITI
            %----------------------------------------------
            HideCursor;
            Screen('FillRect',w,p.foreCol,p.foreRect);
            Screen('FillOval',w,p.fixCol,p.fixRect);
            Screen('DrawingFinished',w);

            % flip ITI
            [time.iti.vblstamp(t) time.iti.onset(t) time.iti.flipstamp(t) time.iti.missed(t)] = Screen('Flip',w);

            if p.portCode == 1
                io64(ioObj,address,p.itiTrigger);
            end

            % Draw fixation in advance
            HideCursor;
            Screen('FillRect',w,p.foreCol,p.foreRect);
            Screen('FillOval',w,p.fixCol,p.fixRect);
            Screen('DrawingFinished', w);

            if p.portCode == 1
                WaitSecs(0.001);
                io64(ioObj,address,p.triggerReset);
            end

            % wait ITI
            WaitSecs('UntilTime',time.iti.onset(t) + stim.ITI(t) - 0.5 * p.refreshCycle);

            %------------------------------------------------------------------
            % fixation
            %------------------------------------------------------------------

            % flip fixation
            [time.fixation.vblstamp(t) time.fixation.onset(t) time.fixation.flipstamp(t) time.fixation.missed(t)] = Screen('Flip',w);

            if p.portCode == 1
                io64(ioObj,address,p.fixTrigger);
            end

            if p.eyeLinkCode == 1
                Eyelink('message', ['Encode' num2str(t+p.nTrials*(b-1))]);
            end

            % Draw sample stimulus in advance
            Screen('FillRect',w,p.foreCol,p.foreRect);
            Screen('FillPoly',w,stim.pColor(t,:), squeeze(stim.triangleOrdin(:,:,t)));
            Screen('FillArc',w,stim.pColor(t,:),stim.arcRectOrdin(t,:),(stim.angle(t,:)+180),180);
            Screen('DrawingFinished',w);

            if p.portCode == 1
                WaitSecs(0.001);
                io64(ioObj,address,p.triggerReset);
            end


            % wait fixation
            WaitSecs('UntilTime',time.fixation.onset(t) + p.fixDur - 0.5 * p.refreshCycle);

            %------------------------------------------------------------------
            % stimuli
            %------------------------------------------------------------------

            % Flip stim
            [time.stim.vblstamp(t) time.stim.onset(t) time.stim.flipstamp(t) time.stim.missed(t)] = Screen('Flip',w);

            if p.portCode == 1
                io64(ioObj,address,p.stimTrigger);
            end

            % Draw delay in advance
            Screen('FillRect',w,p.foreCol,p.foreRect);
            Screen('FillOval',w,p.fixCol,p.fixRect);
            Screen('DrawingFinished',w);

            if p.portCode == 1
                WaitSecs(0.001);
                io64(ioObj,address,p.triggerReset);
            end

            % wait stimuli
            WaitSecs('UntilTime',time.stim.onset(t) + p.stimDur - 0.5 * p.refreshCycle);

            %------------------------------------------------------------------
            % delay
            %------------------------------------------------------------------

            % Flip delay
            [time.delay.vblstamp(t) time.delay.onset(t) time.delay.flipstamp(t) time.delay.missed(t)] = Screen('Flip',w);

            if p.portCode == 1
                io64(ioObj,address,p.delayTrigger);
                WaitSecs(0.001);
                io64(ioObj,address,p.triggerReset);
            end

            WaitSecs('UntilTime',time.delay.onset(t) + p.delay - 0.5 * p.refreshCycle);

            %------------------------------------------------------------------
            % Get response
            %------------------------------------------------------------------

            % Draw ring in advance
            Screen('FillRect',w,p.foreCol,p.foreRect);
            Screen('FillOval',w,p.fixCol,p.fixRect)
            Screen('FrameOval',w,p.rimCol,p.ringPos,p.rimThick,p.rimThick);
            Screen('DrawingFinished',w);

            % Flip ring
            [time.orientationResp.vblstamp(t) time.orientationResp.onset(t) time.orientationResp.flipstamp(t) time.orientationResp.missed(t)] = Screen('Flip',w);

            if p.portCode == 1
                io64(ioObj,address,p.ringTrigger);
            end

            % get response
            SetMouse(p.xCenter+1,p.yCenter+1,w);
            ShowCursor('Arrow');

            if p.portCode == 1
                WaitSecs(0.001);
                io64(ioObj,address,p.triggerReset);
            end

            while GetSecs < time.orientationResp.onset(t) + p.waitOrienResp

                [keyIsDown,secs,keyCode]=KbCheck;

                if keyIsDown
                    kp = find(keyCode);
                    if kp == 27
                        Screen('CloseAll')
                        return
                    end
                end

                [x,y,buttons] = GetMouse(w);

                if any(buttons)
                    rtEnd = GetSecs;
                    if p.portCode == 1
                        io64(ioObj,address,p.orinRepTrigger);
                    end

                    stim.rLoc(t,1)=x; stim.rLoc(t,2)=y;
                    stim.rt(t) = rtEnd-time.orientationResp.onset(t);
                    HideCursor;

                    if p.portCode == 1
                        WaitSecs(0.001);
                        io64(ioObj,address,p.triggerReset);
                    end

                    break;
                end
                if any(buttons)
                    Screen('FillOval',w,p.fixColor,p.fixRect);
                    Screen('Flip',w);
                    break;
                end
            end

            HideCursor;

            %------------------------------------------------------------------
            % Evaluate the orientation response
            %------------------------------------------------------------------
            if ~isnan(stim.rt(t))
                % calculating x&y coords in trig unit
                sin_mouse = (stim.rLoc(t,2) - p.yCenter)/sqrt((stim.rLoc(t,1) - p.xCenter)^2+(stim.rLoc(t,2) - p.yCenter)^2);   % based on the radius given here (tLoc-yCenter), this shows y in trigonometric functions
                cos_mouse = (stim.rLoc(t,1) - p.xCenter)/sqrt((stim.rLoc(t,1) - p.xCenter)^2+(stim.rLoc(t,2) - p.yCenter)^2);   % '' for x

                % trig unit is now converted to degree (0deg starts off from the standard position and goes clockwise)
                if sin_mouse > 0 && cos_mouse > 0                           % 0-90degrees
                    mouse_angle = atan(sin_mouse/cos_mouse);
                elseif sin_mouse > 0 && cos_mouse < 0                       % 90-180degrees
                    mouse_angle = atan(sin_mouse/cos_mouse)+pi;
                elseif sin_mouse < 0 && cos_mouse < 0                       % 180-270degrees
                    mouse_angle = atan(sin_mouse/cos_mouse)+pi;
                elseif sin_mouse < 0 && cos_mouse > 0                       % 270-360degrees
                    mouse_angle = atan(sin_mouse/cos_mouse)+2*pi;
                elseif sin_mouse == 0 && cos_mouse == 1                     % 0degree
                    mouse_angle = 0;
                elseif sin_mouse == 1 && cos_mouse == 0                     % 90degrees
                    mouse_angle = pi/2;
                elseif sin_mouse == 0 && cos_mouse == -1                    % 180 degrees
                    mouse_angle = pi;
                elseif sin_mouse == -1 && cos_mouse == 0                    % 270 degrees
                    mouse_angle = 3/2*pi;
                end;

                % Compute original stim angle
                orig_ang = stim.angle(t,:);                % original degree defined above

                % calculate offset of subject response and actual angle
                tAngle = round(orig_ang); stim.tAngle(t) = tAngle;
                rAngle = round(mouse_angle/pi*180); stim.rAngle(t) = rAngle;
                stim.angleOffset(t) = rAngle-tAngle;

                if stim.angleOffset(t) > 180
                    stim.angleOffset(t) = stim.angleOffset(t)-360;
                elseif stim.angleOffset(t) < -180
                    stim.angleOffset(t) = stim.angleOffset(t)+360;
                end

                if p.practice
                    text = num2str(stim.angleOffset(t));
                    tCenter = [p.xCenter-RectWidth(Screen('TextBounds', w, text))/2 p.yCenter-70];
                    Screen('FillRect',w,p.foreCol,p.foreRect);            % Draw the foreground window
                    Screen('DrawText', w, text, tCenter(1), tCenter(2), p.txtCol);
                    Screen('DrawingFinished', w);
                    Screen('Flip', w);
                    WaitSecs(.5);
                end

            else
                if p.practice
                    text = strcat('Response missed !');
                    tCenter = [p.xCenter-RectWidth(Screen('TextBounds', w, text))/2 p.yCenter-70];
                    Screen('FillRect',w,p.foreCol,p.foreRect);            % Draw the foreground window
                    Screen('DrawText', w, text, tCenter(1), tCenter(2), p.txtCol);
                    Screen('DrawingFinished', w);
                    Screen('Flip', w);
                    WaitSecs(.5);
                end
            end

            Screen('FillRect',w,p.foreCol,p.foreRect);
            Screen('FillOval',w,p.fixCol,p.fixRect);
            Screen('DrawingFinished',w);
            Screen('Flip',w);

            Screen('Close');

            %------------------------------------------------------------------
            % save data file at the end of each trial
            %------------------------------------------------------------------
            save(fName,'p','stim','time');

        end  % end of trial loop

        %----------------------------------------------------------------------
        % End of block feedback
        %----------------------------------------------------------------------

        respisnan = double(isnan(stim.angleOffset));
        respisnonan = find(respisnan==0);

        aveOffset =  round(mean(abs(stim.angleOffset(respisnonan))));
        text1 = ['In that block you were ' num2str(aveOffset) ' degrees off on average.'];
        tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, text1))/2 p.yCenter-120];
        text2 = 'Press the spacebar to continue.';
        tCenter2 = [p.xCenter-RectWidth(Screen('TextBounds', w, text2))/2 p.yCenter-70];
        Screen('FillRect',w,p.foreCol,p.foreRect);            % Draw the foreground window
        Screen('FillOval',w,p.fixCol,p.fixRect);              % Draw the fixation point
        Screen('DrawText', w, text1, tCenter1(1), tCenter1(2), p.txtCol);
        Screen('DrawText', w, text2, tCenter2(1), tCenter2(2), p.txtCol);
        Screen('DrawingFinished', w);
        Screen('Flip', w);
        WaitSecs(.5);

        % Wait for a spacebar press to continue
        while 1
            [keyIsDown,secs,keyCode]=KbCheck;
            if keyIsDown
                kp = find(keyCode);
                if kp == 32                 % spacebar
                    break;
                end
                if kp == 27
                    Screen('CloseAll')
                    return
                end
            end
        end
        WaitSecs(.5);
        %----------------------------------------------------------------------
        if p.eyeLinkCode == 1
            if b ~= p.nBlocks

                EyelinkDoTrackerSetup(eL);

                % do a final check of calibration using driftcorrection
                EyelinkDoDriftCorrection(eL);

                Eyelink('StartRecording');

                Eyelink('message', 'Force Calibrate Finished');	 	 % zero-plot time for EDFVIEW

                errorCheck=Eyelink('checkrecording'); 		         % Check recording status

                if(errorCheck~=0)
                    fprintf('Eyelink checked wrong status.\n');
                    cleanup;                                         % cleanup function
                    Eyelink('ShutDown');
                    Screen('CloseAll');
                end

            end
        end
    end      % end of block loop


    % pack up and go home
    ListenChar(0);
    Screen('CloseAll');
    ShowCursor('Arrow');
    if p.eyeLinkCode == 1
        % save eyelink edf
        Eyelink('StopRecording');
        Eyelink('CloseFile');

        try
            fprintf('Receiving data file ''%s''\n',fileName);
            status=Eyelink('ReceiveFile',tempName ,savePath,1);
            if status > 0
                fprintf('ReceiveFile status %d\n ', status);
            end
            if exist(fileName, 'file')==2
                fprintf('Data file ''%s'' can be found in '' %s\n',fileName, pwd);
            end
        catch
            fprintf('Problem receiving data file ''%s''\n',fileName);
        end
        % cd (savePath);
        save(fullfile(savePath,tempName));
        movefile([savePath,'\',tempName,'.edf'],[savePath,'\',fileName,'.edf']);

        Eyelink('ShutDown');
    end

    %% Sub Functions
    function [struct] = buildTrialStructure(cond,p)

    s = fieldnames(cond);
    struct_orien = [];
    struct_color = [];
    reps = p.nReps;
    for i=1:length(s)
        curVal = getfield(cond,s{i});
        c = [];
        d = [];
        for t = 1:curVal
            c = [c, repmat(t,1,reps)];
            d = [d, repmat(i,1,reps)];
        end
        struct_orien = [struct_orien c];
        struct_color = [struct_color d];
    end
    struct = [struct_orien; struct_color]


    %--------------------------------------------------------------------------
    function z = round2(x,y)
    %ROUND2 rounds number to nearest multiple of arbitrary precision.
    %   Z = ROUND2(X,Y) rounds X to nearest multiple of Y.
    %
    %Example 1: round PI to 2 decimal places
    %   >> round2(pi,0.01)
    %   ans =
    %         3.14
    %
    %Example 2: round PI to 4 decimal places
    %   >> round2(pi,1e-4)
    %   ans =
    %         3.1416
    %
    %Example 3: round PI to 8-bit fraction
    %   >> round2(pi,2^-8)
    %   ans =
    %         3.1406
    %
    %Examples 4-6: round PI to other multiples
    %   >> round2(pi,0.05)
    %   ans =
    %         3.15
    %   >> round2(pi,2)
    %   ans =
    %         4
    %   >> round2(pi,5)
    %   ans =
    %         5
    %
    % See also ROUND.

    % defensive programming
    error(nargchk(2,2,nargin))
    error(nargoutchk(0,1,nargout))
    if numel(y)>1
        error('Y must be scalar')
    end

    %
    z = round(x/y)*y;
    %--------------------------------------------------------------------------
    function pix = deg2pix(deg,vDist,pixSize)
    % Convert degrees of visual angle to pixels for easy specification of
    % stimulus size for PsychToolbox. Returns the size of a stimulus in
    % pixels:
    %
    % INPUTS:
    % deg: desired stim size in degrees of visual angle.
    % vDist: viewing distance.
    % pxSize: pixel size (in same units as viewing distance).
    %
    rad = deg2rad(deg); % convert visual angle from degrees to radians
    sz = vDist*tan(rad); % size of stimulus (in same units as vDistand pxSize)
    pix = round(sz/pixSize); % convert to pixels

    %
    %
