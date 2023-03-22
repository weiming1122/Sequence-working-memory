    function sequence_memory_3items_part02
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
    % defAns = {'',num2str(s)};                                                  % fill in some default answers
    box = inputdlg(prompt, 'Enter Subject Information', 1, defAns);            % actually make GUI
    if length(box) == length(defAns)                                           % simple check for enough input, otherwise bail
        p.subNum = num2str(box{1});  p.rndSeed = str2num(box{2});
        rand('state',p.rndSeed);                                               % seed the generator
    else
        return;
    end

    % Determine if Practice
    prompt = {'Practice (1=yes,0=no)'};
    defAns = {''};                                                             % fill in some default answers
    box = inputdlg(prompt, 'Enter Subject Information', 1, defAns);            % actually make GUI
    if length(box) == length(defAns)                                           % simple check for enough input, otherwise bail
        p.practice = str2num(box{1});                                          % seed the generator
    else
        return;
    end

    % Turn of keyboard echoing
    ListenChar(2);

    %--------------------------------------------------------------------------
    % General Experimental Parameters
    %--------------------------------------------------------------------------
    p.portCode = 0;           % (1: EEG recording; 0: no EEG recording)
    p.eyeLinkCode = 0;        % (1: eyeLink recording; 0: no eyeLink recording)

    p.nBlocks = 4;                                     % number of blocks
    p.nItems = 3                                       % number of items
    p.nColors = 3;                                     % number of colors
    p.nOrientations = 12;                              % number of orientations
    p.nReps = 1;                                       % number of repetitions per color and per orientation for each block
    p.nTrials = p.nColors*p.nOrientations*p.nReps;     % number of trials per block
    p.numProbe = 3                                     % number of probe
    p.waitframes = 1                                   % number of waiting for frame during probe period

    % Dimensions for Stims
    p.vDist = 85;                                     % viewing distance (cm)
    p.px = 0.0248;                                     % pixel size (cm) for monitors: 1920*1080, 21.5 inch
    p.rimSize = deg2pix(6,p.vDist,p.px);               % diameter of rim
    p.rimThick = deg2pix(0.2,p.vDist,p.px);            % thickness of rim
    p.fixSize = deg2pix(0.2,p.vDist,p.px);             % fixation size
    p.triangleH1 = deg2pix(3,p.vDist,p.px);          % height of triangle part of stimulus sample
    p.triangleH2 = deg2pix(0.8,p.vDist,p.px);          % bottom of triangle part of stimulus sample
    p.semicircle = deg2pix(0.8,p.vDist,p.px);          % diameter of semicircle part of stimulus sample
    p.DisProbe = deg2pix(3,p.vDist,p.px);            % distance between probe and center of screen
    p.ProbeRadius = deg2pix(1.5,p.vDist,p.px);         % radius of probe

    % Stimulus orientation for forward encoding model (FEM) analysis
    p.nChans = p.nOrientations;
    p.orienChannels = linspace(15,375-360/p.nChans,p.nChans);        % orientation vector: 15:30:345

    % Timing
    p.refreshRate = 60;
    p.refreshCycle = 1/p.refreshRate;
    p.ITI = [1.0:p.refreshCycle:1.5];                  % ITI ranges between 1000 to 1500 ms
    p.fixDurPre = 0.5;                                 % duration of pre-sample fixation
    p.stimDur = 0.2;                                   % duration of sample stimuli
    p.delay = 1.3;                                     % duration of delay
    p.cueDur = 0.5                                     % duration of cue
    p.fixDurPost = [1.3:p.refreshCycle:1.5]            % duration of post-sample fixation (after cue)
    p.probeDur = 5;                                    % duration of probe
    p.waitFixResp = 5;                                 % duration of waiting for fixation dim test
    p.fixDimLastingTime = 0.2                          % duration of fixation dim
    p.waitOrienResp = 7;                               % duration of waiting for response orientation
    p.nLuminace = p.probeDur * p.refreshRate;          % number element of luminance sequence

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

    % origin rule list each block (half repeat and half mirror)
    p.ruleListOrin = [ones(1,p.nTrials/2),2*ones(1,p.nTrials/2)];              % 1:repeat; 2:mirror

    % orign fixation dim list each block (three quarters no dim and one quarter dim)
    p.fixDimOrin = [zeros(1,3*p.nTrials/4),ones(1,p.nTrials/4)];               % 0:no change; 1:change

    if p.portCode == 1
        % trigger number
        p.itiTrigger = 240;
        p.fixPreTrigger = 241;
        p.stim01Trigger = 242;
        p.delay01Trigger = 243;
        p.stim02Trigger = 244;
        p.delay02Trigger = 245;
        p.stim03Trigger = 246;
        p.delay03Trigger = 247;
        p.ruleTrigger = 248;
        p.fixPostTrigger = 249;
        p.probeTtrigger = 250;
        p.fixTestTrigger = 251;
        p.ringTrigger = 252;
        p.orinRep01Trigger = 253;
        p.orinRep02Trigger = 254;
        p.orinRep03Trigger = 255;
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
    p.foreRect = p.sRect;                       % set the foreground/foreRect to the full screen/p.sRect

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
        fileName = ['sub',num2str(p.subNum,'%02d'),'_part02'];
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

        % Build an output directory & check to make sure it doesn't  already exist

        p.root = pwd;

        if ~exist([p.root, '\Subject Data\'], 'dir')
            mkdir([p.root, '\Subject Data\']);
        end

        % Build an output file and check to make sure that it doesn't exist yet either
        fName = [p.root, '\Subject Data\', 'sub',num2str(p.subNum,'%02d'), '_OrienWM_part02_', num2str(b,'%02d'), '.mat'];

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

        time.fixationPre.vblstamp = nan(1,p.nTrials);
        time.fixationPre.onset = nan(1,p.nTrials);
        time.fixationPre.flipstamp = nan(1,p.nTrials);
        time.fixationPre.missed = nan(1,p.nTrials);

        time.stim01.vblstamp = nan(1,p.nTrials);
        time.stim01.onset = nan(1,p.nTrials);
        time.stim01.flipstamp = nan(1,p.nTrials);
        time.stim01.missed = nan(1,p.nTrials);

        time.stim02.vblstamp = nan(1,p.nTrials);
        time.stim02.onset = nan(1,p.nTrials);
        time.stim02.flipstamp = nan(1,p.nTrials);
        time.stim02.missed = nan(1,p.nTrials);

        time.stim03.vblstamp = nan(1,p.nTrials);
        time.stim03.onset = nan(1,p.nTrials);
        time.stim03.flipstamp = nan(1,p.nTrials);
        time.stim03.missed = nan(1,p.nTrials);

        time.delay01.vblstamp = nan(1,p.nTrials);
        time.delay01.onset = nan(1,p.nTrials);
        time.delay01.flipstamp = nan(1,p.nTrials);
        time.delay01.missed = nan(1,p.nTrials);

        time.delay02.vblstamp = nan(1,p.nTrials);
        time.delay02.onset = nan(1,p.nTrials);
        time.delay02.flipstamp = nan(1,p.nTrials);
        time.delay02.missed = nan(1,p.nTrials);

        time.delay03.vblstamp = nan(1,p.nTrials);
        time.delay03.onset = nan(1,p.nTrials);
        time.delay03.flipstamp = nan(1,p.nTrials);
        time.delay03.missed = nan(1,p.nTrials);

        time.rule.vblstamp = nan(1,p.nTrials);
        time.rule.onset = nan(1,p.nTrials);
        time.rule.flipstamp = nan(1,p.nTrials);
        time.rule.missed = nan(1,p.nTrials);

        time.fixationPost.vblstamp = nan(1,p.nTrials);
        time.fixationPost.onset = nan(1,p.nTrials);
        time.fixationPost.flipstamp = nan(1,p.nTrials);
        time.fixationPost.missed = nan(1,p.nTrials);

        time.probe.vbl = nan(p.nTrials,p.nLuminace+1);

        time.fixationResp.vblstamp = nan(1,p.nTrials);
        time.fixationResp.onset = nan(1,p.nTrials);
        time.fixationResp.flipstamp = nan(1,p.nTrials);
        time.fixationResp.missed = nan(1,p.nTrials);

        time.orientationResp.vblstamp = nan(1,p.nTrials);
        time.orientationResp.onset = nan(1,p.nTrials);
        time.orientationResp.flipstamp = nan(1,p.nTrials);
        time.orientationResp.missed = nan(1,p.nTrials);

        %----------------------------------------------------------------------
        % Control parameters
        %----------------------------------------------------------------------

        % preallocation stimuli label vectors
        stim.stim01.orienChanLab = nan(1,p.nTrials);               % stimulus orientation label
        stim.stim01.angle = nan(p.nTrials,1);                      % stimulus orientation angle
        stim.stim01.Degree = nan(p.nTrials,1);                     % stimulus orientation degree
        stim.stim01.A1Ordin =  nan(p.nTrials,2);                   % ordination of first vertex of triangle (part of stimulus)
        stim.stim01.A2Ordin =  nan(p.nTrials,2);                   % ordination of second vertex of triangle (part of stimulus)
        stim.stim01.A3Ordin =  nan(p.nTrials,2);                   % ordination of third vertex of triangle (part of stimulus)
        stim.stim01.triangleOrdin = nan(3,2,p.nTrials);            % ordination of triangle (part of stimulus)
        stim.stim01.arcRectOrdin = nan(p.nTrials,4);               % ordination of rect of semicircle (part of stimulus)
        stim.stim01.color = nan(1,p.nTrials);                      % stimulus color label
        stim.stim01.pColor = nan(p.nTrials,3);                     % stimulus color

        stim.stim02.orienChanLab = nan(1,p.nTrials);               % stimulus orientation label
        stim.stim02.angle = nan(p.nTrials,1);                      % stimulus orientation angle
        stim.stim02.Degree = nan(p.nTrials,1);                     % stimulus orientation degree
        stim.stim02.A1Ordin =  nan(p.nTrials,2);                   % ordination of first vertex of triangle (part of stimulus)
        stim.stim02.A2Ordin =  nan(p.nTrials,2);                   % ordination of second vertex of triangle (part of stimulus)
        stim.stim02.A3Ordin =  nan(p.nTrials,2);                   % ordination of third vertex of triangle (part of stimulus)
        stim.stim02.triangleOrdin = nan(3,2,p.nTrials);            % ordination of triangle (part of stimulus)
        stim.stim02.arcRectOrdin = nan(p.nTrials,4);               % ordination of rect of semicircle (part of stimulus)
        stim.stim02.color = nan(1,p.nTrials);                      % stimulus color label
        stim.stim02.pColor = nan(p.nTrials,3);                     % stimulus color

        stim.stim03.orienChanLab = nan(1,p.nTrials);               % stimulus orientation label
        stim.stim03.angle = nan(p.nTrials,1);                      % stimulus orientation angle
        stim.stim03.Degree = nan(p.nTrials,1);                     % stimulus orientation degree
        stim.stim03.A1Ordin =  nan(p.nTrials,2);                   % ordination of first vertex of triangle (part of stimulus)
        stim.stim03.A2Ordin =  nan(p.nTrials,2);                   % ordination of second vertex of triangle (part of stimulus)
        stim.stim03.A3Ordin =  nan(p.nTrials,2);                   % ordination of third vertex of triangle (part of stimulus)
        stim.stim03.triangleOrdin = nan(3,2,p.nTrials);            % ordination of triangle (part of stimulus)
        stim.stim03.arcRectOrdin = nan(p.nTrials,4);               % ordination of rect of semicircle (part of stimulus)
        stim.stim03.color = nan(1,p.nTrials);                      % stimulus color label
        stim.stim03.pColor = nan(p.nTrials,3);                     % stimulus color

        % preallocation fixation dim response vectors
        stim.fixDim.tResp = nan(1,p.nTrials);           % target value (same as stim.fixDimList: 1 means dim; 0 means no dim)
        stim.fixDim.rResp = nan(1,p.nTrials);           % actual reported value
        stim.fixDim.rt = nan(1,p.nTrials);              % response time

        % preallocation oreintation response vectors
        stim.stim01.tAngle = nan(1,p.nTrials);          % target value (same as stim.probedPos)
        stim.stim01.rAngle = nan(1,p.nTrials);          % actual reported value
        stim.stim01.rt = nan(1,p.nTrials);              % response time
        stim.stim01.rLoc = nan(p.nTrials,2);            % coordinates of mouseclick
        stim.stim01.angleOffset = nan(1,p.nTrials);     % offset error

        stim.stim02.tAngle = nan(1,p.nTrials);          % target value (same as stim.probedPos)
        stim.stim02.rAngle = nan(1,p.nTrials);          % actual reported value
        stim.stim02.rt = nan(1,p.nTrials);              % response time
        stim.stim02.rLoc = nan(p.nTrials,2);            % coordinates of mouseclick
        stim.stim02.angleOffset = nan(1,p.nTrials);     % offset error

        stim.stim03.tAngle = nan(1,p.nTrials);          % target value (same as stim.probedPos)
        stim.stim03.rAngle = nan(1,p.nTrials);          % actual reported value
        stim.stim03.rt = nan(1,p.nTrials);              % response time
        stim.stim03.rLoc = nan(p.nTrials,2);            % coordinates of mouseclick
        stim.stim03.angleOffset = nan(1,p.nTrials);     % offset error

        % preallocation probe location
        stim.probe.locAngStar = nan(p.nTrials,1);             % start angle of first probe
        stim.probe.locAngEnd = nan(p.nTrials,1);              % end angle of first probe
        stim.probe.locAngVect = nan(p.nTrials,4);             % angle vector of three probes
        stim.probe.locDegVect = nan(p.nTrials,4);             % degree vector of three probes
        stim.probe.xPosVector = nan(p.nTrials,4);             % x coordination vector of three probes
        stim.probe.yPosVector = nan(p.nTrials,4);             % y coordination vector of three probes
        stim.probe.locProbe01 = nan(p.nTrials,4);             % coordination of rect of first probe
        stim.probe.locProbe02 = nan(p.nTrials,4);             % coordination of rect of second probe
        stim.probe.locProbe03 = nan(p.nTrials,4);             % coordination of rect of third probe
        stim.probe.locProbeAll = nan(3,4,p.nTrials);          % coordination of rect of three probes
        stim.probe.posRandID = nan(p.nTrials,3);              % suffle the positions of three probes
        stim.probe.locProbeAllRand = nan(3,4,p.nTrials);      % coordination of rect of three probes after suffle

        % preallocation probe luminance
        stim.probe.lumSeq = nan(3,p.nLuminace,p.nTrials);                % all item luminance sequence
        stim.stim01.colorlumSeq = nan(p.nLuminace,3,p.nTrials);          % item1 color luminance sequence
        stim.stim02.colorlumSeq = nan(p.nLuminace,3,p.nTrials);          % item2 color luminance sequence
        stim.stim03.colorlumSeq = nan(p.nLuminace,3,p.nTrials);          % item3 color luminance sequence

        % create orientation Channel sequence and color sequence for stimuli

        [tStructure] = buildTrialStructure(p);

        stim.stim01.orienChanLab = mod(tStructure(1,:),12);                    % 1:1:12 (15:30:345)
        stim.stim01.orienChanLab(find(stim.stim01.orienChanLab==0))=12;        % change 0 to 12
        stim.stim01.color = ceil(tStructure(1,:)/12);                          % 1:red; 2:green; 3:blue

        stim.stim02.orienChanLab = mod(tStructure(2,:),12);                    % 1:1:12 (15:30:345)
        stim.stim02.orienChanLab(find(stim.stim02.orienChanLab==0))=12;        % change 0 to 12
        stim.stim02.color = ceil(tStructure(2,:)/12);                          % 1:red; 2:green; 3:blue

        stim.stim03.orienChanLab = mod(tStructure(3,:),12);                    % 1:1:12 (15:30:345)
        stim.stim03.orienChanLab(find(stim.stim03.orienChanLab==0))=12;        % change 0 to 12
        stim.stim03.color = ceil(tStructure(3,:)/12);                          % 1:red; 2:green; 3:blue

        % create rule list per block
        stim.ruleList = p.ruleListOrin(randperm(p.nTrials));

        % create fixation change list per block
        stim.fixDimList = p.fixDimOrin(randperm(p.nTrials));
        stim.fixDim.tResp = stim.fixDimList;

        %----------------------------------------------------------------------
        % Define Stimulus Parameters Before starting block
        %----------------------------------------------------------------------

        cnt=1;

        for t = 1:p.nTrials

            if rem(cnt,p.nTrials) == 0

                % Display for Block loading
                text1 = ['Block ' num2str(b,'%02d') ' Loading'];
                tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, text1))/2 p.yCenter-120];
                text2 = [num2str(round2(cnt/p.nTrials,1e-2)*100), ' %'];
                tCenter2 = [p.xCenter-RectWidth(Screen('TextBounds', w, text2))/2 p.yCenter-70];

                Screen('FillRect',w,p.foreCol,p.foreRect);                             % Draw the foreground window
                Screen('FillOval',w,p.fixCol,p.fixRect);                               % Draw the fixation point
                Screen('DrawText', w, text1, tCenter1(1), tCenter1(2), p.txtCol);
                Screen('DrawText', w, text2, tCenter2(1), tCenter2(2), p.txtCol);
                Screen('DrawingFinished', w);
                Screen('Flip', w);

                [keyIsDown,secs,keyCode]=KbCheck;

                if keyIsDown
                    kp = find(keyCode); kp=kp(1);
                    if kp == 27               % 'Esc'
                        Screen('CloseAll')
                        return
                    end
                end
            end

            % set ITI for trial
            stim.ITI(t) = randsample(p.ITI,1);

            % set duration of post sample fixation
            stim.fixDurPost(t) = randsample(p.fixDurPost,1);

            % set fixation dim start time
            stim.fixDimStart(t) = rand(1)*2.6 + 1;

            % define location of stimuli
            stim.stim01.angle(t,:) = p.orienChannels(stim.stim01.orienChanLab(t));
            stim.stim01.Degree(t,:) = stim.stim01.angle(t,:)*pi/180;
            stim.stim01.A1Ordin(t,:) = [(p.xCenter + p.triangleH1 * cos(stim.stim01.Degree(t,:))),(p.yCenter + p.triangleH1 * sin(stim.stim01.Degree(t,:)))];
            stim.stim01.A2Ordin(t,:) = [(p.xCenter + p.triangleH2 * cos(stim.stim01.Degree(t,:) + 90*pi/180)),(p.yCenter + p.triangleH2 * sin(stim.stim01.Degree(t,:) + 90*pi/180))];
            stim.stim01.A3Ordin(t,:) = [(p.xCenter + p.triangleH2 * cos(stim.stim01.Degree(t,:) - 90*pi/180)),(p.yCenter + p.triangleH2 * sin(stim.stim01.Degree(t,:) - 90*pi/180))];
            stim.stim01.triangleOrdin(:,:,t) = [ stim.stim01.A1Ordin(t,:); stim.stim01.A2Ordin(t,:); stim.stim01.A3Ordin(t,:)];
            stim.stim01.arcRectOrdin(t,:) = [(p.xCenter - p.semicircle), (p.yCenter - p.semicircle), (p.xCenter + p.semicircle),(p.yCenter + p.semicircle)];

            stim.stim02.angle(t,:) = p.orienChannels(stim.stim02.orienChanLab(t));
            stim.stim02.Degree(t,:) = stim.stim02.angle(t,:)*pi/180;
            stim.stim02.A1Ordin(t,:) = [(p.xCenter + p.triangleH1 * cos(stim.stim02.Degree(t,:))),(p.yCenter + p.triangleH1 * sin(stim.stim02.Degree(t,:)))];
            stim.stim02.A2Ordin(t,:) = [(p.xCenter + p.triangleH2 * cos(stim.stim02.Degree(t,:) + 90*pi/180)),(p.yCenter + p.triangleH2 * sin(stim.stim02.Degree(t,:) + 90*pi/180))];
            stim.stim02.A3Ordin(t,:) = [(p.xCenter + p.triangleH2 * cos(stim.stim02.Degree(t,:) - 90*pi/180)),(p.yCenter + p.triangleH2 * sin(stim.stim02.Degree(t,:) - 90*pi/180))];
            stim.stim02.triangleOrdin(:,:,t) = [ stim.stim02.A1Ordin(t,:); stim.stim02.A2Ordin(t,:); stim.stim02.A3Ordin(t,:)];
            stim.stim02.arcRectOrdin(t,:) = [(p.xCenter - p.semicircle), (p.yCenter - p.semicircle), (p.xCenter + p.semicircle),(p.yCenter + p.semicircle)];

            stim.stim03.angle(t,:) = p.orienChannels(stim.stim03.orienChanLab(t));
            stim.stim03.Degree(t,:) = stim.stim03.angle(t,:)*pi/180;
            stim.stim03.A1Ordin(t,:) = [(p.xCenter + p.triangleH1 * cos(stim.stim03.Degree(t,:))),(p.yCenter + p.triangleH1 * sin(stim.stim03.Degree(t,:)))];
            stim.stim03.A2Ordin(t,:) = [(p.xCenter + p.triangleH2 * cos(stim.stim03.Degree(t,:) + 90*pi/180)),(p.yCenter + p.triangleH2 * sin(stim.stim03.Degree(t,:) + 90*pi/180))];
            stim.stim03.A3Ordin(t,:) = [(p.xCenter + p.triangleH2 * cos(stim.stim03.Degree(t,:) - 90*pi/180)),(p.yCenter + p.triangleH2 * sin(stim.stim03.Degree(t,:) - 90*pi/180))];
            stim.stim03.triangleOrdin(:,:,t) = [ stim.stim03.A1Ordin(t,:); stim.stim03.A2Ordin(t,:); stim.stim03.A3Ordin(t,:)];
            stim.stim03.arcRectOrdin(t,:) = [(p.xCenter - p.semicircle), (p.yCenter - p.semicircle), (p.xCenter + p.semicircle),(p.yCenter + p.semicircle)];

            % define colors of stimuli
            if stim.stim01.color(t) == 1
                stim.stim01.pColor(t,:) = p.red;    % 1:red
            elseif stim.stim01.color(t) == 2
                stim.stim01.pColor(t,:) = p.green;  % 2:green
            else
                stim.stim01.pColor(t,:) = p.blue;   % 3:blue
            end

            if stim.stim02.color(t) == 1
                stim.stim02.pColor(t,:) = p.red;    % 1:red
            elseif stim.stim02.color(t) == 2
                stim.stim02.pColor(t,:) = p.green;  % 2:green
            else
                stim.stim02.pColor(t,:) = p.blue;   % 3:blue
            end

            if stim.stim03.color(t) == 1
                stim.stim03.pColor(t,:) = p.red;    % 1:red
            elseif stim.stim03.color(t) == 2
                stim.stim03.pColor(t,:) = p.green;  % 2:green
            else
                stim.stim03.pColor(t,:) = p.blue;   % 3:blue
            end

            % set probe location
            stim.probe.locAngStar(t,:) = randi(121);
            stim.probe.locAngEnd(t,:) = stim.probe.locAngStar(t,:) + 360;
            stim.probe.locAngVect(t,:) = linspace(stim.probe.locAngStar(t,:), stim.probe.locAngEnd(t,:), p.numProbe+1);
            stim.probe.locDegVect(t,:) = stim.probe.locAngVect(t,:) * pi/180;
            stim.probe.xPosVector(t,:) = cos(stim.probe.locDegVect(t,:)) .* p.DisProbe + p.xCenter;
            stim.probe.yPosVector(t,:) = sin(stim.probe.locDegVect(t,:)) .* p.DisProbe + p.yCenter;

            stim.probe.locProbe01(t,:) = CenterRectOnPointd([0,0,2*p.ProbeRadius,2*p.ProbeRadius], stim.probe.xPosVector(t,1), stim.probe.yPosVector(t,1));
            stim.probe.locProbe02(t,:) = CenterRectOnPointd([0,0,2*p.ProbeRadius,2*p.ProbeRadius], stim.probe.xPosVector(t,2), stim.probe.yPosVector(t,2));
            stim.probe.locProbe03(t,:) = CenterRectOnPointd([0,0,2*p.ProbeRadius,2*p.ProbeRadius], stim.probe.xPosVector(t,3), stim.probe.yPosVector(t,3));

            stim.probe.locProbeAll(:,:,t)=[stim.probe.locProbe01(t,:);stim.probe.locProbe02(t,:);stim.probe.locProbe03(t,:)];
            stim.probe.posRandID(t,:)= randperm(3);
            stim.probe.locProbeAllRand(:,:,t)= stim.probe.locProbeAll(stim.probe.posRandID(t,:),:,t);

            % creat three independent luminance sequence for three stimi colors
            for k = 1:3
                lumr=(rand([1,p.nLuminace])*1);     % 300 luminance sequence
                alum=lumr-mean(lumr);               % normalization
                N=length(alum);
                y=fft(alum,N);
                y(find(real(y)==0))=1+0i;
                % nl=0:N-1; % align to matlab index
                % f=nl*100/N;% sequence of frequency
                mag=abs(y);
                realp=real(y);
                imagp=imag(y);
                weig=sqrt((repmat((mean(mag))^2,1,N))./(realp.^2+imagp.^2));   % compute weight to normalize to mean. the case of mag=0 will lead error
                unistim=(weig.*realp)+(weig.*imagp)*1i;
                ulummid=ifft(unistim);
                ulummid= ulummid+mean(lumr);
                stim.probe.lumSeq(k,:,t) = ulummid;                            % first row for item1; second row for item2; third row for item3
            end

            stim.stim01.colorlumSeq(:,:,t) = (stim.stim01.pColor(t,:)'*stim.probe.lumSeq(1,:,t))';  % item1 color luminance sequence
            stim.stim02.colorlumSeq(:,:,t) = (stim.stim02.pColor(t,:)'*stim.probe.lumSeq(2,:,t))';  % item2 color luminance sequence
            stim.stim03.colorlumSeq(:,:,t) = (stim.stim03.pColor(t,:)'*stim.probe.lumSeq(3,:,t))';  % item3 color luminance sequence

            cnt=cnt+1;
        end

        % Display for Block loading
        text1 = ['You have ' num2str(p.nBlocks-(b-1)) ' blocks remaining.'];
        tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, text1))/2 p.yCenter-120];
        text2 = 'Press the spacebar to continue';
        tCenter2 = [p.xCenter-RectWidth(Screen('TextBounds', w, text2))/2 p.yCenter-70];
        Screen('FillRect',w,p.foreCol,p.foreRect);              % Draw the foreground window
        Screen('FillOval',w,p.fixCol,p.fixRect);                % Draw the fixation point
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

            %------------------------------------------------------------------
            % ITI
            %------------------------------------------------------------------
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
            % fixation pre
            %------------------------------------------------------------------

            % flip fixation pre
            [time.fixationPre.vblstamp(t) time.fixationPre.onset(t) time.fixationPre.flipstamp(t) time.fixationPre.missed(t)] = Screen('Flip',w);

            if p.portCode == 1
                io64(ioObj,address,p.fixPreTrigger);
            end

            if p.eyeLinkCode == 1
                Eyelink('message', ['Encode' num2str(t+p.nTrials*(b-1))]);
            end

            % Draw sample stimulus01 in advance
            Screen('FillRect',w,p.foreCol,p.foreRect);
            Screen('FillPoly',w,stim.stim01.pColor(t,:), squeeze(stim.stim01.triangleOrdin(:,:,t)));
            Screen('FillArc',w,stim.stim01.pColor(t,:),stim.stim01.arcRectOrdin(t,:),(stim.stim01.angle(t,:)+180),180);
            Screen('DrawingFinished',w);

            if p.portCode == 1
                WaitSecs(0.001);
                io64(ioObj,address,p.triggerReset);
            end

            % wait fixation
            WaitSecs('UntilTime',time.fixationPre.onset(t) + p.fixDurPre - 0.5 * p.refreshCycle);

            %------------------------------------------------------------------
            % stimuli 01
            %------------------------------------------------------------------

            % Flip stim 01
            [time.stim01.vblstamp(t) time.stim01.onset(t) time.stim01.flipstamp(t) time.stim01.missed(t)] = Screen('Flip',w);

            if p.portCode == 1
                io64(ioObj,address,p.stim01Trigger);
            end

            % Draw delay01 in advance
            Screen('FillRect',w,p.foreCol,p.foreRect);
            Screen('FillOval',w,p.fixCol,p.fixRect);
            Screen('DrawingFinished',w);

            if p.portCode == 1
                WaitSecs(0.001);
                io64(ioObj,address,p.triggerReset);
            end

            % wait stimuli01
            WaitSecs('UntilTime',time.stim01.onset(t) + p.stimDur - 0.5 * p.refreshCycle);

            %------------------------------------------------------------------
            % delay01
            %------------------------------------------------------------------

            % Flip delay01
            [time.delay01.vblstamp(t) time.delay01.onset(t) time.delay01.flipstamp(t) time.delay01.missed(t)] = Screen('Flip',w);

            if p.portCode == 1
                io64(ioObj,address,p.delay01Trigger);
            end

            % Draw sample stimulus02 in advance
            Screen('FillRect',w,p.foreCol,p.foreRect);
            Screen('FillPoly',w,stim.stim02.pColor(t,:), squeeze(stim.stim02.triangleOrdin(:,:,t)));
            Screen('FillArc',w,stim.stim02.pColor(t,:),stim.stim02.arcRectOrdin(t,:),(stim.stim02.angle(t,:)+180),180);
            Screen('DrawingFinished',w);

            if p.portCode == 1
                WaitSecs(0.001);
                io64(ioObj,address,p.triggerReset);
            end

            % wait delay01
            WaitSecs('UntilTime',time.delay01.onset(t) + p.delay - 0.5 * p.refreshCycle);

            %------------------------------------------------------------------
            % stimuli 02
            %------------------------------------------------------------------

            % Flip stim 02
            [time.stim02.vblstamp(t) time.stim02.onset(t) time.stim02.flipstamp(t) time.stim02.missed(t)] = Screen('Flip',w);

            if p.portCode == 1
                io64(ioObj,address,p.stim02Trigger);
            end

            % Draw delay02 in advance
            Screen('FillRect',w,p.foreCol,p.foreRect);
            Screen('FillOval',w,p.fixCol,p.fixRect);
            Screen('DrawingFinished',w);

            if p.portCode == 1
                WaitSecs(0.001);
                io64(ioObj,address,p.triggerReset);
            end

            % wait stimuli02
            WaitSecs('UntilTime',time.stim02.onset(t) + p.stimDur - 0.5 * p.refreshCycle);

            %------------------------------------------------------------------
            % delay02
            %------------------------------------------------------------------

            % Flip delay02
            [time.delay02.vblstamp(t) time.delay02.onset(t) time.delay02.flipstamp(t) time.delay02.missed(t)] = Screen('Flip',w);

            if p.portCode == 1
                io64(ioObj,address,p.delay02Trigger);
            end

            % Draw sample stimulus03 in advance
            Screen('FillRect',w,p.foreCol,p.foreRect);
            Screen('FillPoly',w,stim.stim03.pColor(t,:), squeeze(stim.stim03.triangleOrdin(:,:,t)));
            Screen('FillArc',w,stim.stim03.pColor(t,:),stim.stim03.arcRectOrdin(t,:),(stim.stim03.angle(t,:)+180),180);
            Screen('DrawingFinished',w);

            if p.portCode == 1
                WaitSecs(0.001);
                io64(ioObj,address,p.triggerReset);
            end

            % wait delay02
            WaitSecs('UntilTime',time.delay02.onset(t) + p.delay-0.5 * p.refreshCycle);

            %------------------------------------------------------------------
            % stimuli 03
            %------------------------------------------------------------------

            % Flip stim 03
            [time.stim03.vblstamp(t) time.stim03.onset(t) time.stim03.flipstamp(t) time.stim03.missed(t)] = Screen('Flip',w);

            if p.portCode == 1
                io64(ioObj,address,p.stim03Trigger);
            end

            % Draw delay03 in advance
            Screen('FillRect',w,p.foreCol,p.foreRect);
            Screen('FillOval',w,p.fixCol,p.fixRect);
            Screen('DrawingFinished',w);

            if p.portCode == 1
                WaitSecs(0.001);
                io64(ioObj,address,p.triggerReset);
            end

            % wait stimuli03
            WaitSecs('UntilTime',time.stim03.onset(t) + p.stimDur - 0.5 * p.refreshCycle);

            %------------------------------------------------------------------
            % delay03
            %------------------------------------------------------------------

            % Flip delay03
            [time.delay03.vblstamp(t) time.delay03.onset(t) time.delay03.flipstamp(t) time.delay03.missed(t)] = Screen('Flip',w);

            if p.portCode == 1
                io64(ioObj,address,p.delay03Trigger);
            end

            % Draw Rule Cue in advance
            if stim.ruleList(t) == 1
                text = 'repeat';
                tCenter = [p.xCenter-RectWidth(Screen('TextBounds', w, text))/2 p.yCenter-20];
                Screen('FillRect',w,p.foreCol,p.foreRect);                     % Draw the foreground window
                Screen('DrawText', w, text, tCenter(1), tCenter(2), p.txtCol);
                Screen('DrawingFinished', w);
            else
                text = 'Mirror';
                tCenter = [p.xCenter-RectWidth(Screen('TextBounds', w, text))/2 p.yCenter-20];
                Screen('FillRect',w,p.foreCol,p.foreRect);                     % Draw the foreground window
                Screen('DrawText', w, text, tCenter(1), tCenter(2), p.txtCol);
                Screen('DrawingFinished', w);
            end

            if p.portCode == 1
                WaitSecs(0.001);
                io64(ioObj,address,p.triggerReset);
            end

            % wait delay03
            WaitSecs('UntilTime',time.delay03.onset(t) + p.delay - 0.5 * p.refreshCycle);

            %------------------------------------------------------------------
            % Rule Cue
            %------------------------------------------------------------------

            % Flip Rule Cue
            [time.rule.vblstamp(t) time.rule.onset(t) time.rule.flipstamp(t) time.rule.missed(t)] = Screen('Flip',w);

            if p.portCode == 1
                io64(ioObj,address,p.ruleTrigger);
            end

            % Draw Fixtion in advance
            Screen('FillRect',w,p.foreCol,p.foreRect);
            Screen('FillOval',w,p.fixCol,p.fixRect);
            Screen('DrawingFinished', w);

            if p.portCode == 1
                WaitSecs(0.001);
                io64(ioObj,address,p.triggerReset);
            end

            % wait Rule Cue
            WaitSecs('UntilTime',time.rule.onset(t) + p.cueDur - 0.5 * p.refreshCycle);

            %------------------------------------------------------------------
            % fixation post
            %------------------------------------------------------------------

            % flip fixation post
            [time.fixationPost.vblstamp(t) time.fixationPost.onset(t) time.fixationPost.flipstamp(t) time.fixationPost.missed(t)] = Screen('Flip',w);

            if p.portCode == 1
                io64(ioObj,address,p.fixPostTrigger);
                WaitSecs(0.001);
                io64(ioObj,address,p.triggerReset);
            end

            % wait fixation post
            WaitSecs('UntilTime',time.fixationPost.onset(t) + stim.fixDurPost(t) - 0.5*p.refreshCycle);

            %------------------------------------------------------------------
            % probe
            %------------------------------------------------------------------

            time.probe.vbl(t,1) = GetSecs;
            i = 1;

            if p.portCode == 1
                io64(ioObj,address,p.probeTtrigger);
            end
            if p.eyeLinkCode == 1
                Eyelink('message', ['Probe' num2str(t+p.nTrials*(b-1))]);
            end
            while time.probe.vbl(t,i) < time.probe.vbl(t,1) + p.probeDur && i <= 300

                Screen('FillOval', w, stim.stim01.colorlumSeq(i,:,t), stim.probe.locProbeAllRand(1,:,t));
                Screen('FillOval', w, stim.stim02.colorlumSeq(i,:,t), stim.probe.locProbeAllRand(2,:,t));
                Screen('FillOval', w, stim.stim03.colorlumSeq(i,:,t), stim.probe.locProbeAllRand(3,:,t));

                % fix dim or not
                Screen('FillOval',w,p.fixCol,p.fixRect);
                if stim.fixDimList(t) == 1
                    if GetSecs > time.probe.vbl(t,1) + stim.fixDimStart(t) && GetSecs < time.probe.vbl(t,1) + stim.fixDimStart(t) + p.fixDimLastingTime
                        Screen('FillOval',w,0.6*p.fixCol,p.fixRect);
                    end
                end

                i = i+1;
                time.probe.vbl(t,i) = Screen('Flip', w, time.probe.vbl(t,i-1) + (p.waitframes - 0.5) * p.refreshCycle);
            end

            if p.portCode == 1
                io64(ioObj,address,p.triggerReset);
            end

            %------------------------------------------------------------------
            % fix test
            %------------------------------------------------------------------

            text = '?';
            tCenter = [p.xCenter-RectWidth(Screen('TextBounds', w, text))/2 p.yCenter - 20];
            Screen('FillRect',w,p.foreCol,p.foreRect);                     % Draw the foreground window
            Screen('DrawText', w, text, tCenter(1), tCenter(2), p.txtCol);
            Screen('DrawingFinished', w);


            [time.fixationResp.vblstamp(t) time.fixationResp.onset(t) time.fixationResp.flipstamp(t) time.fixationResp.missed(t)] = Screen('Flip',w);

            if p.portCode == 1
                io64(ioObj,address,p.fixTestTrigger);
                WaitSecs(0.001);
                io64(ioObj,address,p.triggerReset);
            end

            while GetSecs < time.fixationResp.onset(t) + p.waitFixResp

                [keyIsDown, secs, keyCode] = KbCheck;
                [x, y, buttonsFix] = GetMouse(w);

                if buttonsFix(1)
                    fprintf('left click 1\n');
                    stim.fixDim.rResp(1, t) = 1;
                    stim.fixDim.rt(1, t) = GetSecs - time.fixationResp.onset(t);
                    break
                elseif buttonsFix(3)
                    fprintf('right click 1\n');
                    stim.fixDim.rResp(1, t) = 0;
                    stim.fixDim.rt(1, t) = GetSecs - time.fixationResp.onset(t);
                    break
                elseif keyIsDown
                    kp = find(keyCode);
                    if kp == 27
                        Screen('CloseAll')
                        return
                    end
                end
            end

            while any(buttonsFix) % if already down, wait for release
                [x, y, buttonsFix] = GetMouse(w);
            end

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
            clicks = 0;

            if p.portCode == 1
                WaitSecs(0.001);
                io64(ioObj,address,p.triggerReset);
            end

            while GetSecs < time.orientationResp.onset(t) + p.waitOrienResp

                [keyIsDown,secs,keyCode] = KbCheck;

                if keyIsDown
                    kp = find(keyCode);
                    if kp == 27
                        Screen('CloseAll')
                        return
                    end
                end

                [x,y,buttonsOrien] = GetMouse(w);

                if sum(buttonsOrien) > 0
                    clicks = clicks + 1;
                    if sum(buttonsOrien) > 0 && clicks == 1
                        rtEnd01 = GetSecs;
                        if p.portCode == 1
                            io64(ioObj,address,p.orinRep01Trigger);
                        end
                        stim.stim01.rt(t) = rtEnd01-time.orientationResp.onset(t);
                        stim.stim01.rLoc(t,1)=x; stim.stim01.rLoc(t,2)=y;
                        SetMouse(p.xCenter+1,p.yCenter+1,w);
                        if p.portCode == 1
                            WaitSecs(0.001);
                            io64(ioObj,address,p.triggerReset);
                        end
                    end

                    while sum(buttonsOrien) > 0 && clicks == 1
                        [x,y,buttonsOrien] = GetMouse(w);
                    end

                    if sum(buttonsOrien) > 0 && clicks == 2
                        rtEnd02 = GetSecs;
                        if p.portCode == 1
                            io64(ioObj,address,p.orinRep02Trigger);
                        end
                        stim.stim02.rt(t) = rtEnd02-rtEnd01;
                        stim.stim02.rLoc(t,1)=x; stim.stim02.rLoc(t,2)=y;
                        SetMouse(p.xCenter+1,p.yCenter+1,w);
                        if p.portCode == 1
                            WaitSecs(0.001);
                            io64(ioObj,address,p.triggerReset);
                        end
                    end

                    while sum(buttonsOrien) > 0 && clicks == 2
                        [x,y,buttonsOrien] = GetMouse(w);
                    end

                    if sum(buttonsOrien) > 0 && clicks == 3
                        rtEnd03 = GetSecs;
                        if p.portCode == 1
                            io64(ioObj,address,p.orinRep03Trigger);
                        end
                        stim.stim03.rt(t) = rtEnd03-rtEnd02;
                        stim.stim03.rLoc(t,1)=x; stim.stim03.rLoc(t,2)=y;
                        if p.portCode == 1
                            WaitSecs(0.001);
                            io64(ioObj,address,p.triggerReset);
                        end
                        break;
                    end
                end
            end

            HideCursor;

            %------------------------------------------------------------------
            % Evaluate the orientation response
            %------------------------------------------------------------------
            if ~isnan(stim.stim01.rt(t)) && ~isnan(stim.stim02.rt(t)) && ~isnan(stim.stim03.rt(t))      % as the case response to all three orientations
                % calculating x&y coords in trig unit
                sin_mouse01 = (stim.stim01.rLoc(t,2) - p.yCenter)/sqrt((stim.stim01.rLoc(t,1) - p.xCenter)^2+(stim.stim01.rLoc(t,2) - p.yCenter)^2);   % based on the radius given here (tLoc-yCenter), this shows y in trigonometric functions
                cos_mouse01 = (stim.stim01.rLoc(t,1) - p.xCenter)/sqrt((stim.stim01.rLoc(t,1) - p.xCenter)^2+(stim.stim01.rLoc(t,2) - p.yCenter)^2);   % '' for x
                sin_mouse02 = (stim.stim02.rLoc(t,2) - p.yCenter)/sqrt((stim.stim02.rLoc(t,1) - p.xCenter)^2+(stim.stim02.rLoc(t,2) - p.yCenter)^2);   % based on the radius given here (tLoc-yCenter), this shows y in trigonometric functions
                cos_mouse02 = (stim.stim02.rLoc(t,1) - p.xCenter)/sqrt((stim.stim02.rLoc(t,1) - p.xCenter)^2+(stim.stim02.rLoc(t,2) - p.yCenter)^2);   % '' for x
                sin_mouse03 = (stim.stim03.rLoc(t,2) - p.yCenter)/sqrt((stim.stim03.rLoc(t,1) - p.xCenter)^2+(stim.stim03.rLoc(t,2) - p.yCenter)^2);   % based on the radius given here (tLoc-yCenter), this shows y in trigonometric functions
                cos_mouse03 = (stim.stim03.rLoc(t,1) - p.xCenter)/sqrt((stim.stim03.rLoc(t,1) - p.xCenter)^2+(stim.stim03.rLoc(t,2) - p.yCenter)^2);   % '' for x

                % trig unit is now converted to degree (0deg starts off from the standard position and goes clockwise)
                if sin_mouse01 > 0 && cos_mouse01 > 0                          % 0-90degrees
                    mouse_angle01 = atan(sin_mouse01/cos_mouse01);
                elseif sin_mouse01 > 0 && cos_mouse01 < 0                      % 90-180degrees
                    mouse_angle01 = atan(sin_mouse01/cos_mouse01)+pi;
                elseif sin_mouse01 < 0 && cos_mouse01 < 0                      % 180-270degrees
                    mouse_angle01 = atan(sin_mouse01/cos_mouse01)+pi;
                elseif sin_mouse01 < 0 && cos_mouse01 > 0                      % 270-360degrees
                    mouse_angle01 = atan(sin_mouse01/cos_mouse01)+2*pi;
                elseif sin_mouse01 == 0 && cos_mouse01 == 1                    % 0degree
                    mouse_angle01 = 0;
                elseif sin_mouse01 == 1 && cos_mouse01 == 0                    % 90degrees
                    mouse_angle01 = pi/2;
                elseif sin_mouse01 == 0 && cos_mouse01 == -1                   % 180 degrees
                    mouse_angle01 = pi;
                elseif sin_mouse01 == -1 && cos_mouse01 == 0                   % 270 degrees
                    mouse_angle01 = 3/2*pi;
                end;

                if sin_mouse02 > 0 && cos_mouse02 > 0                          % 0-90degrees
                    mouse_angle02 = atan(sin_mouse02/cos_mouse02);
                elseif sin_mouse02 > 0 && cos_mouse02 < 0                      % 90-180degrees
                    mouse_angle02 = atan(sin_mouse02/cos_mouse02)+pi;
                elseif sin_mouse02 < 0 && cos_mouse02 < 0                      % 180-270degrees
                    mouse_angle02 = atan(sin_mouse02/cos_mouse02)+pi;
                elseif sin_mouse02 < 0 && cos_mouse02 > 0                      % 270-360degrees
                    mouse_angle02 = atan(sin_mouse02/cos_mouse02)+2*pi;
                elseif sin_mouse02 == 0 && cos_mouse02 == 1                    % 0degree
                    mouse_angle02 = 0;
                elseif sin_mouse02 == 1 && cos_mouse02 == 0                    % 90degrees
                    mouse_angle02 = pi/2;
                elseif sin_mouse02 == 0 && cos_mouse02 == -1                   % 180 degrees
                    mouse_angle02 = pi;
                elseif sin_mouse02 == -1 && cos_mouse02 == 0                   % 270 degrees
                    mouse_angle02 = 3/2*pi;
                end;

                if sin_mouse03 > 0 && cos_mouse03 > 0                          % 0-90degrees
                    mouse_angle03 = atan(sin_mouse03/cos_mouse03);
                elseif sin_mouse03 > 0 && cos_mouse03 < 0                      % 90-180degrees
                    mouse_angle03 = atan(sin_mouse03/cos_mouse03)+pi;
                elseif sin_mouse03 < 0 && cos_mouse03 < 0                      % 180-270degrees
                    mouse_angle03 = atan(sin_mouse03/cos_mouse03)+pi;
                elseif sin_mouse03 < 0 && cos_mouse03 > 0                      % 270-360degrees
                    mouse_angle03 = atan(sin_mouse03/cos_mouse03)+2*pi;
                elseif sin_mouse03 == 0 && cos_mouse03 == 1                    % 0degree
                    mouse_angle03 = 0;
                elseif sin_mouse03 == 1 && cos_mouse03 == 0                    % 90degrees
                    mouse_angle03 = pi/2;
                elseif sin_mouse03 == 0 && cos_mouse03 == -1                    % 180 degrees
                    mouse_angle03 = pi;
                elseif sin_mouse03 == -1 && cos_mouse03 == 0                    % 270 degrees
                    mouse_angle03 = 3/2*pi;
                end;

                % Compute original stim angle
                orig_ang01 = stim.stim01.angle(t,:);                % original degree defined above
                orig_ang02 = stim.stim02.angle(t,:);                % original degree defined above
                orig_ang03 = stim.stim03.angle(t,:);                % original degree defined above

                if stim.ruleList(t) == 1                   % repeat trial
                    % calculate offset of subject response and actual angle
                    tAngle01 = round(orig_ang01); stim.stim01.tAngle(t) = tAngle01;
                    rAngle01 = round(mouse_angle01/pi*180); stim.stim01.rAngle(t) = rAngle01;
                    stim.stim01.angleOffset(t) = rAngle01-tAngle01;

                    tAngle02 = round(orig_ang02); stim.stim02.tAngle(t) = tAngle02;
                    rAngle02 = round(mouse_angle02/pi*180); stim.stim02.rAngle(t) = rAngle02;
                    stim.stim02.angleOffset(t) = rAngle02-tAngle02;

                    tAngle03 = round(orig_ang03); stim.stim03.tAngle(t) = tAngle03;
                    rAngle03 = round(mouse_angle03/pi*180); stim.stim03.rAngle(t) = rAngle03;
                    stim.stim03.angleOffset(t) = rAngle03-tAngle03;
                else                                        % mirror trial
                    % calculate offset of subject response and actual angle
                    tAngle01 = round(orig_ang03); stim.stim01.tAngle(t) = tAngle01;
                    rAngle01 = round(mouse_angle01/pi*180); stim.stim01.rAngle(t) = rAngle01;
                    stim.stim01.angleOffset(t) = rAngle01-tAngle01;

                    tAngle02 = round(orig_ang02); stim.stim02.tAngle(t) = tAngle02;
                    rAngle02 = round(mouse_angle02/pi*180); stim.stim02.rAngle(t) = rAngle02;
                    stim.stim02.angleOffset(t) = rAngle02-tAngle02;

                    tAngle03 = round(orig_ang01); stim.stim03.tAngle(t) = tAngle03;
                    rAngle03 = round(mouse_angle03/pi*180); stim.stim03.rAngle(t) = rAngle03;
                    stim.stim03.angleOffset(t) = rAngle03-tAngle03;
                end

                if stim.stim01.angleOffset(t) > 180
                    stim.stim01.angleOffset(t) = stim.stim01.angleOffset(t)-360;
                elseif stim.stim01.angleOffset(t) < -180
                    stim.stim01.angleOffset(t) = stim.stim01.angleOffset(t)+360;
                end

                if stim.stim02.angleOffset(t) > 180
                    stim.stim02.angleOffset(t) = stim.stim02.angleOffset(t)-360;
                elseif stim.stim02.angleOffset(t) < -180
                    stim.stim02.angleOffset(t) = stim.stim02.angleOffset(t)+360;
                end

                if stim.stim03.angleOffset(t) > 180
                    stim.stim03.angleOffset(t) = stim.stim03.angleOffset(t)-360;
                elseif stim.stim03.angleOffset(t) < -180
                    stim.stim03.angleOffset(t) = stim.stim03.angleOffset(t)+360;
                end

                if p.practice
                    text = strcat(num2str(stim.stim01.angleOffset(t)),';',num2str(stim.stim02.angleOffset(t)),';',num2str(stim.stim03.angleOffset(t)));
                    tCenter = [p.xCenter-RectWidth(Screen('TextBounds', w, text))/2 p.yCenter-70];
                    Screen('FillRect',w,p.foreCol,p.foreRect);            % Draw the foreground window
                    Screen('DrawText', w, text, tCenter(1), tCenter(2), p.txtCol);
                    Screen('DrawingFinished', w);
                    Screen('Flip', w);
                    WaitSecs(.5);
                end

            else                                                 % as the case not response to all three orientations
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
        isnan01 = double(isnan(stim.stim01.angleOffset));
        isnan02 = double(isnan(stim.stim02.angleOffset));
        isnan03 = double(isnan(stim.stim03.angleOffset));
        isnonan01 = find(isnan01==0);
        isnonan02 = find(isnan02==0);
        isnonan03 = find(isnan03==0);
        isnonan = intersect(intersect(isnonan01,isnonan02),isnonan03);

        aveOffset =  round(mean(abs(stim.stim01.angleOffset(isnonan))+abs(stim.stim02.angleOffset(isnonan))+abs(stim.stim03.angleOffset(isnonan))));
        text1 = ['In that block you were ' num2str(aveOffset) ' degrees off on average across three items.'];
        tCenter1 = [p.xCenter-RectWidth(Screen('TextBounds', w, text1))/2 p.yCenter-120];
        text2 = 'Press the spacebar to continue.';
        tCenter2 = [p.xCenter-RectWidth(Screen('TextBounds', w, text2))/2 p.yCenter-70];
        Screen('FillRect',w,p.foreCol,p.foreRect);              % Draw the foreground window
        Screen('FillOval',w,p.fixCol,p.fixRect);                % Draw the fixation point
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
                if kp == 32        % spacebar
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

    end          % end of block loop


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
    function [struct] = buildTrialStructure(p)
    while 1
        seq = [1:1:p.nTrials];
        S01 = seq;S02 = seq;S03 = seq;
        Stim01 = S01(randperm(p.nTrials));
        for i = 1: p.nTrials
            if ((Stim01(i)<=12) & (Stim01(i)>0))
                seq_out02 = [1:1:12,Stim01(i)+12,Stim01(i)+24];
                seq_remain02 = setdiff(S02,seq_out02);
                if isempty(seq_remain02)
                    break
                end
                Stim02(i) = seq_remain02(randi(length(seq_remain02)));
                S02 = setdiff(S02,Stim02(i));
                if ((Stim02(i)<=24) & (Stim02(i)>12))
                    seq_out03 = union(seq_out02,[13:1:24,Stim02(i)+12]);
                    seq_remain03 = setdiff(S03,seq_out03);
                    if isempty(seq_remain03)
                        break
                    end
                    Stim03(i) = seq_remain03(randi(length(seq_remain03)));
                    S03 = setdiff(S03,Stim03(i));
                end
                if ((Stim02(i)<=36) & (Stim02(i)>24))
                    seq_out03 = union(seq_out02,[25:1:36,Stim02(i)-12]);
                    seq_remain03 = setdiff(S03,seq_out03);
                    if isempty(seq_remain03)
                        break
                    end
                    Stim03(i) = seq_remain03(randi(length(seq_remain03)));
                    S03 = setdiff(S03,Stim03(i));
                end
            end

            if ((Stim01(i)<=24) & (Stim01(i)>12))
                seq_out02 = [13:1:24,Stim01(i)-12,Stim01(i)+12];
                seq_remain02 = setdiff(S02,seq_out02);
                if isempty(seq_remain02)
                    break
                end
                Stim02(i) = seq_remain02(randi(length(seq_remain02)));
                S02 = setdiff(S02,Stim02(i));
                if ((Stim02(i)<=12) & (Stim02(i)>0))
                    seq_out03 = union(seq_out02,[1:1:12,Stim02(i)+24]);
                    seq_remain03 = setdiff(S03,seq_out03);
                    if isempty(seq_remain03)
                        break
                    end
                    Stim03(i) = seq_remain03(randi(length(seq_remain03)));
                    S03 = setdiff(S03,Stim03(i));
                end
                if ((Stim02(i)<=36) & (Stim02(i)>24))
                    seq_out03 = union(seq_out02,[25:1:36,Stim02(i)-24]);
                    seq_remain03 = setdiff(S03,seq_out03);
                    if isempty(seq_remain03)
                        break
                    end
                    Stim03(i) = seq_remain03(randi(length(seq_remain03)));
                    S03 = setdiff(S03,Stim03(i));
                end
            end
            if ((Stim01(i)<=36) & (Stim01(i)>24))
                seq_out02 = [25:1:36,Stim01(i)-12,Stim01(i)-24];
                seq_remain02 = setdiff(S02,seq_out02);
                if isempty(seq_remain02)
                    break
                end
                Stim02(i) = seq_remain02(randi(length(seq_remain02)));
                S02 = setdiff(S02,Stim02(i));
                if ((Stim02(i)<=12) & (Stim02(i)>0))
                    seq_out03 = union(seq_out02,[1:1:12,Stim02(i)+12]);
                    seq_remain03 = setdiff(S03,seq_out03);
                    if isempty(seq_remain03)
                        break
                    end
                    Stim03(i) = seq_remain03(randi(length(seq_remain03)));
                    S03 = setdiff(S03,Stim03(i));
                end
                if ((Stim02(i)<=24) & (Stim02(i)>12))
                    seq_out03 = union(seq_out02,[13:1:24,Stim02(i)-12]);
                    seq_remain03 = setdiff(S03,seq_out03);
                    if isempty(seq_remain03)
                        break
                    end
                    Stim03(i) = seq_remain03(randi(length(seq_remain03)));
                    S03 = setdiff(S03,Stim03(i));
                end
            end
        end
        if ~isempty(seq_remain02)&&~isempty(seq_remain03)
            break
        end
    end
    struct = [Stim01;Stim02;Stim03];

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
