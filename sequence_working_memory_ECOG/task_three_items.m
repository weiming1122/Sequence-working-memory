function sequence_memory_3items_clinical
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

p.nBlocks = 4;                                    % number of blocks
p.nItems = 3;                                     % number of items
p.nOrientations = 8;                              % number of orientations
p.nReps = 3;                                      % number of repetitions per orientation for each block
p.nTrials = p.nOrientations * p.nReps;

% Dimensions for Stims
p.vDist = 85;                                      % viewing distance (cm)
p.px = 0.0248;                                     % pixel size (cm) for monitors: 1920*1080, 21.5 inch
p.rimSize = deg2pix(6,p.vDist,p.px);               % diameter of rim
p.rimThick = deg2pix(0.2,p.vDist,p.px);            % thickness of rim
p.fixSize = deg2pix(0.2,p.vDist,p.px);             % fixation size
p.triangleH1 = deg2pix(3,p.vDist,p.px);            % height of triangle part of stimulus sample
p.triangleH2 = deg2pix(0.8,p.vDist,p.px);          % bottom of triangle part of stimulus sample
p.semicircle = deg2pix(0.8,p.vDist,p.px);          % diameter of semicircle part of stimulus sample


% Stimulus orientation for forward encoding model (FEM) analysis
p.nChans = p.nOrientations;
p.orienChannels = linspace(0,360-360/p.nChans,p.nChans);        % orientation vector: 0:45:315

% Timing
p.refreshRate = 60;
p.refreshCycle = 1/p.refreshRate;
p.ITI = [1.5:p.refreshCycle:2.0];                  % ITI ranges between 1500 to 2000 ms
p.fixDurPre = 0.5;                                 % duration of pre-sample fixation
p.stimDur = 0.5;                                   % duration of sample stimuli
p.delay = 1.5;                                     % duration of delay
p.cueDur = 1.0;                                     % duration of cue
p.fixDurPost = 1.5;                                 % duration of post-sample fixation (after cue)
p.waitOrienResp = 10;                               % duration of waiting for response orientation

% Color information
p.black = [0 0 0];                     % color of black
p.dimgray = [90 90 90];                % color of Dimgray
p.gray = [128 128 128];                % color of gray

p.txtCol = p.black
p.foreCol = p.gray;
p.fixCol = p.dimgray;
p.stimColor = p.dimgray;
p.rimCol = p.dimgray;

% origin rule list each block (half repeat and half mirror)
p.ruleListOrin = [ones(1,p.nTrials/2),2*ones(1,p.nTrials/2)];              % 1:repeat; 2:mirror

if p.portCode == 1
    % trigger number
    p.stim01Trigger = 10;
    p.ruleTrigger = 20;
    p.orinRep01Trigger = 30;
    p.triggerReset = 0;
end

%--------------------------------------------------------------------------
% Build the stimulus display
%--------------------------------------------------------------------------
Screen('Preference','VisualDebugLevel', 3);
AssertOpenGL;                                                   % Make sure we're on an OpenGL-compatible machine
s = max(Screen('Screens'));                                     % Grab a handle for the display ID (should return 0 for single-monitor setups)
Screen('Preference', 'SkipSyncTests', 0);
% [w,p.sRect] = Screen('OpenWindow',s,p.foreCol);                     % full-screen
% [w,p.sRect] = Screen('OpenWindow',s,p.foreCol,[1,0,1599,900]);      % when connect with other monitor
[w,p.sRect] = Screen('OpenWindow',s,p.foreCol,[1,0,1919,1080]);
p.ifi = Screen('GetFlipInterval', w);                           % estimate the monitor flip interval
HideCursor;
Priority(MaxPriority(w));                                       % set priority to max to discourage interruptions
Screen(w,'TextSize',[40]);
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
    ioObj = io32; % initialize inpoutx32.dll device driver
    status = io32(ioObj);
    if status == 0
        disp('inpoutx32.dll successfully installed.')
    else
        error('inpoutx32.dll installation failed.')
    end
    address = hex2dec('C100');   % LPT3 Status port addr
    io32(ioObj,address,0);
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
    
    time.rim.vblstamp = nan(1,p.nTrials);
    time.rim.onset = nan(1,p.nTrials);
    time.rim.flipstamp = nan(1,p.nTrials);
    time.rim.missed = nan(1,p.nTrials);
    
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
    
    stim.stim02.orienChanLab = nan(1,p.nTrials);               % stimulus orientation label
    stim.stim02.angle = nan(p.nTrials,1);                      % stimulus orientation angle
    stim.stim02.Degree = nan(p.nTrials,1);                     % stimulus orientation degree
    stim.stim02.A1Ordin =  nan(p.nTrials,2);                   % ordination of first vertex of triangle (part of stimulus)
    stim.stim02.A2Ordin =  nan(p.nTrials,2);                   % ordination of second vertex of triangle (part of stimulus)
    stim.stim02.A3Ordin =  nan(p.nTrials,2);                   % ordination of third vertex of triangle (part of stimulus)
    stim.stim02.triangleOrdin = nan(3,2,p.nTrials);            % ordination of triangle (part of stimulus)
    stim.stim02.arcRectOrdin = nan(p.nTrials,4);               % ordination of rect of semicircle (part of stimulus)
    
    stim.stim03.orienChanLab = nan(1,p.nTrials);               % stimulus orientation label
    stim.stim03.angle = nan(p.nTrials,1);                      % stimulus orientation angle
    stim.stim03.Degree = nan(p.nTrials,1);                     % stimulus orientation degree
    stim.stim03.A1Ordin =  nan(p.nTrials,2);                   % ordination of first vertex of triangle (part of stimulus)
    stim.stim03.A2Ordin =  nan(p.nTrials,2);                   % ordination of second vertex of triangle (part of stimulus)
    stim.stim03.A3Ordin =  nan(p.nTrials,2);                   % ordination of third vertex of triangle (part of stimulus)
    stim.stim03.triangleOrdin = nan(3,2,p.nTrials);            % ordination of triangle (part of stimulus)
    stim.stim03.arcRectOrdin = nan(p.nTrials,4);               % ordination of rect of semicircle (part of stimulus)
    
    % preallocation oreintation response vectors
    stim.resp01.tAngle = nan(1,p.nTrials);          % target value (same as stim.probedPos)
    stim.resp01.rAngle = nan(1,p.nTrials);          % actual reported value
    stim.resp01.rt = nan(1,p.nTrials);              % response time
    stim.resp01.rLoc = nan(p.nTrials,2);            % coordinates of mouseclick
    stim.resp01.angleOffset = nan(1,p.nTrials);     % offset error
    
    stim.resp02.tAngle = nan(1,p.nTrials);          % target value (same as stim.probedPos)
    stim.resp02.rAngle = nan(1,p.nTrials);          % actual reported value
    stim.resp02.rt = nan(1,p.nTrials);              % response time
    stim.resp02.rLoc = nan(p.nTrials,2);            % coordinates of mouseclick
    stim.resp02.angleOffset = nan(1,p.nTrials);     % offset error
    
    stim.resp03.tAngle = nan(1,p.nTrials);          % target value (same as stim.probedPos)
    stim.resp03.rAngle = nan(1,p.nTrials);          % actual reported value
    stim.resp03.rt = nan(1,p.nTrials);              % response time
    stim.resp03.rLoc = nan(p.nTrials,2);            % coordinates of mouseclick
    stim.resp03.angleOffset = nan(1,p.nTrials);     % offset error
    
    % create orientation Channel sequence and color sequence for stimuli
    
    [tStructure] = buildTrialStructure(p);
    
    stim.stim01.orienChanLab = tStructure(1,:)                     % 1:1:8 (0:45:315)
    stim.stim02.orienChanLab = tStructure(2,:);                    % 1:1:8 (0:45:315)
    stim.stim03.orienChanLab = tStructure(3,:);                    % 1:1:8 (0:45:315)
    
    % create rule list per block
    stim.ruleList = p.ruleListOrin(randperm(p.nTrials));
    
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
        
        
        % Draw fixation in advance
        HideCursor;
        Screen('FillRect',w,p.foreCol,p.foreRect);
        Screen('FillOval',w,p.fixCol,p.fixRect);
        Screen('DrawingFinished', w);
        
        %------------------------------------------------------------------
        % fixation pre
        %------------------------------------------------------------------
        
        % flip fixation pre
        [time.fixationPre.vblstamp(t) time.fixationPre.onset(t) time.fixationPre.flipstamp(t) time.fixationPre.missed(t)] = ... 
            Screen('Flip',w,time.iti.onset(t) + stim.ITI(t) - 1 * p.refreshCycle);
        
        % Draw sample stimulus01 in advance
        Screen('FillRect',w,p.foreCol,p.foreRect);
        Screen('FillPoly',w,p.stimColor, squeeze(stim.stim01.triangleOrdin(:,:,t)));
        Screen('FillArc',w,p.stimColor,stim.stim01.arcRectOrdin(t,:),(stim.stim01.angle(t,:)+180),180);
        Screen('DrawingFinished',w);
        
        %------------------------------------------------------------------
        % stimuli 01
        %------------------------------------------------------------------
        
        % Flip stim 01
        [time.stim01.vblstamp(t) time.stim01.onset(t) time.stim01.flipstamp(t) time.stim01.missed(t)] = ... 
            Screen('Flip',w,time.fixationPre.onset(t) + p.fixDurPre - 1 * p.refreshCycle);
        
        if p.portCode == 1
            io32(ioObj,address,p.stim01Trigger);
        end
        
        % Draw delay01 in advance
        Screen('FillRect',w,p.foreCol,p.foreRect);
        Screen('FillOval',w,p.fixCol,p.fixRect);
        Screen('DrawingFinished',w);
        
        if p.portCode == 1
            WaitSecs(0.1);
            io32(ioObj,address,p.triggerReset);
        end
        
        %------------------------------------------------------------------
        % delay01
        %------------------------------------------------------------------
        
        % Flip delay01
        [time.delay01.vblstamp(t) time.delay01.onset(t) time.delay01.flipstamp(t) time.delay01.missed(t)] = ... 
            Screen('Flip',w,time.stim01.onset(t)+ p.stimDur- 1 * p.refreshCycle);
        
        % Draw sample stimulus02 in advance
        Screen('FillRect',w,p.foreCol,p.foreRect);
        Screen('FillPoly',w,p.stimColor, squeeze(stim.stim02.triangleOrdin(:,:,t)));
        Screen('FillArc',w,p.stimColor,stim.stim02.arcRectOrdin(t,:),(stim.stim02.angle(t,:)+180),180);
        Screen('DrawingFinished',w);
        
        %------------------------------------------------------------------
        % stimuli 02
        %------------------------------------------------------------------
        
        % Flip stim 02
        [time.stim02.vblstamp(t) time.stim02.onset(t) time.stim02.flipstamp(t) time.stim02.missed(t)] = ... 
            Screen('Flip',w,time.delay01.onset(t) + p.delay - 1 * p.refreshCycle);
        
        % Draw delay02 in advance
        Screen('FillRect',w,p.foreCol,p.foreRect);
        Screen('FillOval',w,p.fixCol,p.fixRect);
        Screen('DrawingFinished',w);
        
        %------------------------------------------------------------------
        % delay02
        %------------------------------------------------------------------
        
        % Flip delay02
        [time.delay02.vblstamp(t) time.delay02.onset(t) time.delay02.flipstamp(t) time.delay02.missed(t)] = ...
            Screen('Flip',w,time.stim02.onset(t) + p.stimDur- 1 * p.refreshCycle);
        
        % Draw sample stimulus03 in advance
        Screen('FillRect',w,p.foreCol,p.foreRect);
        Screen('FillPoly',w,p.stimColor, squeeze(stim.stim03.triangleOrdin(:,:,t)));
        Screen('FillArc',w,p.stimColor,stim.stim03.arcRectOrdin(t,:),(stim.stim03.angle(t,:)+180),180);
        Screen('DrawingFinished',w);
        
        %------------------------------------------------------------------
        % stimuli 03
        %------------------------------------------------------------------
        
        % Flip stim 03
        [time.stim03.vblstamp(t) time.stim03.onset(t) time.stim03.flipstamp(t) time.stim03.missed(t)] = ...
            Screen('Flip',w,time.delay02.onset(t) + p.delay- 1 * p.refreshCycle);
        
        % Draw delay03 in advance
        Screen('FillRect',w,p.foreCol,p.foreRect);
        Screen('FillOval',w,p.fixCol,p.fixRect);
        Screen('DrawingFinished',w);
        
        %------------------------------------------------------------------
        % delay03
        %------------------------------------------------------------------
        
        % Flip delay03
        [time.delay03.vblstamp(t) time.delay03.onset(t) time.delay03.flipstamp(t) time.delay03.missed(t)] = ...
            Screen('Flip',w,time.stim03.onset(t) + p.stimDur - 1 * p.refreshCycle);
        
        % Draw Rule Cue in advance
        if stim.ruleList(t) == 1
            text = 'repeat';
%             text = double('ÕýÐò');
            tCenter = [p.xCenter-RectWidth(Screen('TextBounds', w, text))/2 p.yCenter-20];
            Screen('FillRect',w,p.foreCol,p.foreRect);                     % Draw the foreground window
            Screen('DrawText', w, text, tCenter(1), tCenter(2), p.txtCol);
            Screen('DrawingFinished', w);
        else
             text = 'Mirror';
%             text = double('µ¹Ðò');
            tCenter = [p.xCenter-RectWidth(Screen('TextBounds', w, text))/2 p.yCenter-20];
            Screen('FillRect',w,p.foreCol,p.foreRect);                     % Draw the foreground window
            Screen('DrawText', w, text, tCenter(1), tCenter(2), p.txtCol);
            Screen('DrawingFinished', w);
        end
        
        %------------------------------------------------------------------
        % Rule Cue
        %------------------------------------------------------------------
        
        % Flip Rule Cue
        [time.rule.vblstamp(t) time.rule.onset(t) time.rule.flipstamp(t) time.rule.missed(t)] = ...
            Screen('Flip',w,time.delay03.onset(t) + p.delay - 1 * p.refreshCycle);
        
        if p.portCode == 1
            io32(ioObj,address,p.ruleTrigger);
        end
        
        % Draw Fixtion in advance
        Screen('FillRect',w,p.foreCol,p.foreRect);
        Screen('FillOval',w,p.fixCol,p.fixRect);
        Screen('DrawingFinished', w);
        
        if p.portCode == 1
            WaitSecs(0.1);
            io32(ioObj,address,p.triggerReset);
        end
        
        %------------------------------------------------------------------
        % fixation post
        %------------------------------------------------------------------
        
        % flip fixation post
        [time.fixationPost.vblstamp(t) time.fixationPost.onset(t) time.fixationPost.flipstamp(t) time.fixationPost.missed(t)] = ...
            Screen('Flip',w,time.rule.onset(t) + p.cueDur - 1 * p.refreshCycle);
        
        % Draw ring in advance
        Screen('FillRect',w,p.foreCol,p.foreRect);
        Screen('FillOval',w,p.fixCol,p.fixRect)
        Screen('FrameOval',w,p.rimCol,p.ringPos,p.rimThick,p.rimThick);
        Screen('DrawingFinished',w);
        
        %------------------------------------------------------------------
        % Get response
        %------------------------------------------------------------------
        
        % Flip ring
        [time.rim.vblstamp(t) time.rim.onset(t) time.rim.flipstamp(t) time.rim.missed(t)] = ... 
            Screen('Flip',w,time.fixationPost.onset(t) + p.fixDurPost - 1*p.refreshCycle);
        
        % get response
        SetMouse(p.xCenter+1,p.yCenter+1,w);
        % SetMouse(799+1,430+1,w);
        ShowCursor('Arrow');
        clicks = 0;
        
        while GetSecs < time.rim.onset(t) + p.waitOrienResp
            
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
                        io32(ioObj,address,p.orinRep01Trigger);
                    end
                    stim.resp01.rt(t) = rtEnd01-time.rim.onset(t);
                    stim.resp01.rLoc(t,1)=x; stim.resp01.rLoc(t,2)=y;
                    SetMouse(p.xCenter+1,p.yCenter+1,w);
                    if p.portCode == 1
                        WaitSecs(0.1);
                        io32(ioObj,address,p.triggerReset);
                    end
                end
                
                while sum(buttonsOrien) > 0 && clicks == 1
                    [x,y,buttonsOrien] = GetMouse(w);
                end
                
                if sum(buttonsOrien) > 0 && clicks == 2
                    rtEnd02 = GetSecs;
                    stim.resp02.rt(t) = rtEnd02-rtEnd01;
                    stim.resp02.rLoc(t,1)=x; stim.resp02.rLoc(t,2)=y;
                    SetMouse(p.xCenter+1,p.yCenter+1,w);
                end
                
                while sum(buttonsOrien) > 0 && clicks == 2
                    [x,y,buttonsOrien] = GetMouse(w);
                end
                
                if sum(buttonsOrien) > 0 && clicks == 3
                    rtEnd03 = GetSecs;
                    stim.resp03.rt(t) = rtEnd03-rtEnd02;
                    stim.resp03.rLoc(t,1)=x; stim.resp03.rLoc(t,2)=y;
                    break;
                end
            end
        end
        
        HideCursor;
        
        %------------------------------------------------------------------
        % Evaluate the orientation response
        %------------------------------------------------------------------
        if ~isnan(stim.resp01.rt(t)) && ~isnan(stim.resp02.rt(t)) && ~isnan(stim.resp03.rt(t))      % as the case response to all three orientations
            % calculating x&y coords in trig unit
            sin_mouse01 = (stim.resp01.rLoc(t,2) - p.yCenter)/sqrt((stim.resp01.rLoc(t,1) - p.xCenter)^2+(stim.resp01.rLoc(t,2) - p.yCenter)^2);   % based on the radius given here (tLoc-yCenter), this shows y in trigonometric functions
            cos_mouse01 = (stim.resp01.rLoc(t,1) - p.xCenter)/sqrt((stim.resp01.rLoc(t,1) - p.xCenter)^2+(stim.resp01.rLoc(t,2) - p.yCenter)^2);   % '' for x
            sin_mouse02 = (stim.resp02.rLoc(t,2) - p.yCenter)/sqrt((stim.resp02.rLoc(t,1) - p.xCenter)^2+(stim.resp02.rLoc(t,2) - p.yCenter)^2);   % based on the radius given here (tLoc-yCenter), this shows y in trigonometric functions
            cos_mouse02 = (stim.resp02.rLoc(t,1) - p.xCenter)/sqrt((stim.resp02.rLoc(t,1) - p.xCenter)^2+(stim.resp02.rLoc(t,2) - p.yCenter)^2);   % '' for x
            sin_mouse03 = (stim.resp03.rLoc(t,2) - p.yCenter)/sqrt((stim.resp03.rLoc(t,1) - p.xCenter)^2+(stim.resp03.rLoc(t,2) - p.yCenter)^2);   % based on the radius given here (tLoc-yCenter), this shows y in trigonometric functions
            cos_mouse03 = (stim.resp03.rLoc(t,1) - p.xCenter)/sqrt((stim.resp03.rLoc(t,1) - p.xCenter)^2+(stim.resp03.rLoc(t,2) - p.yCenter)^2);   % '' for x
            
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
                tAngle01 = round(orig_ang01); stim.resp01.tAngle(t) = tAngle01;
                rAngle01 = round(mouse_angle01/pi*180); stim.resp01.rAngle(t) = rAngle01;
                stim.resp01.angleOffset(t) = rAngle01-tAngle01;
                
                tAngle02 = round(orig_ang02); stim.resp02.tAngle(t) = tAngle02;
                rAngle02 = round(mouse_angle02/pi*180); stim.resp02.rAngle(t) = rAngle02;
                stim.resp02.angleOffset(t) = rAngle02-tAngle02;
                
                tAngle03 = round(orig_ang03); stim.resp03.tAngle(t) = tAngle03;
                rAngle03 = round(mouse_angle03/pi*180); stim.resp03.rAngle(t) = rAngle03;
                stim.resp03.angleOffset(t) = rAngle03-tAngle03;
            else                                        % mirror trial
                % calculate offset of subject response and actual angle
                tAngle01 = round(orig_ang03); stim.resp01.tAngle(t) = tAngle01;
                rAngle01 = round(mouse_angle01/pi*180); stim.resp01.rAngle(t) = rAngle01;
                stim.resp01.angleOffset(t) = rAngle01-tAngle01;
                
                tAngle02 = round(orig_ang02); stim.resp02.tAngle(t) = tAngle02;
                rAngle02 = round(mouse_angle02/pi*180); stim.resp02.rAngle(t) = rAngle02;
                stim.resp02.angleOffset(t) = rAngle02-tAngle02;
                
                tAngle03 = round(orig_ang01); stim.resp03.tAngle(t) = tAngle03;
                rAngle03 = round(mouse_angle03/pi*180); stim.resp03.rAngle(t) = rAngle03;
                stim.resp03.angleOffset(t) = rAngle03-tAngle03;
            end
            
            if stim.resp01.angleOffset(t) > 180
                stim.resp01.angleOffset(t) = stim.resp01.angleOffset(t)-360;
            elseif stim.resp01.angleOffset(t) < -180
                stim.resp01.angleOffset(t) = stim.resp01.angleOffset(t)+360;
            end
            
            if stim.resp02.angleOffset(t) > 180
                stim.resp02.angleOffset(t) = stim.resp02.angleOffset(t)-360;
            elseif stim.resp02.angleOffset(t) < -180
                stim.resp02.angleOffset(t) = stim.resp02.angleOffset(t)+360;
            end
            
            if stim.resp03.angleOffset(t) > 180
                stim.resp03.angleOffset(t) = stim.resp03.angleOffset(t)-360;
            elseif stim.resp03.angleOffset(t) < -180
                stim.resp03.angleOffset(t) = stim.resp03.angleOffset(t)+360;
            end
            
            if p.practice
                text = strcat(num2str(stim.resp01.angleOffset(t)),';  ',13,num2str(stim.resp02.angleOffset(t)),';  ',13, ...
                    num2str(stim.resp03.angleOffset(t)));
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
    isnan01 = double(isnan(stim.resp01.angleOffset));
    isnan02 = double(isnan(stim.resp02.angleOffset));
    isnan03 = double(isnan(stim.resp03.angleOffset));
    isnonan01 = find(isnan01==0);
    isnonan02 = find(isnan02==0);
    isnonan03 = find(isnan03==0);
    isnonan = intersect(intersect(isnonan01,isnonan02),isnonan03);
    
    aveOffset =  round(mean(abs(stim.resp01.angleOffset(isnonan))+abs(stim.resp02.angleOffset(isnonan))+abs(stim.resp03.angleOffset(isnonan))));
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
    
end          % end of block loop


% pack up and go home
ListenChar(0);
Screen('CloseAll');
ShowCursor('Arrow');


%% Sub Functions
function [struct] = buildTrialStructure(p)
struct = [];
for r= 1: p.nReps
    while 1
        seq = [1:1:p.nOrientations];
        S01 = seq;S02 = seq;S03 = seq;
        Stim01 = S01(randperm(p.nOrientations));
        for i = 1: p.nOrientations
            seq_out02 = Stim01(i);
            seq_remain02 = setdiff(S02,seq_out02);
            if isempty(seq_remain02)
                break;
            end
            Stim02(i) = seq_remain02(randi(length(seq_remain02)));
            S02 = setdiff(S02,Stim02(i));
            seq_out03 = union(seq_out02,Stim02(i));
            seq_remain03 = setdiff(S03,seq_out03);
            if isempty(seq_remain03)
                break;
            end
            Stim03(i) = seq_remain03(randi(length(seq_remain03)));
            S03 = setdiff(S03,Stim03(i));
        end
        if ~isempty(seq_remain02)&&~isempty(seq_remain03)
            break
        end
    end
    struct_temp = [Stim01;Stim02;Stim03];
    struct = cat(2,struct,struct_temp);
end

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
