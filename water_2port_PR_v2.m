 %{
----------------------------------------------------------------------------

This file is part of the Sanworks Bpod repository
Copyright (C) 2016 Sanworks LLC, Sound Beach, New York, USA

----------------------------------------------------------------------------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3.

This program is distributed  WITHOUT ANY WARRANTY and without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}
function Sarah
% This protocol is a starting point for a task with three ports.
% Here animals just learn that water is available at three adjacent ports.
% Water is availabe when lights are on.
% Written by Sarah Starosta and Aubrey Siebels 10/2018.
%
% SETUP
% You will need:
% - A Bpod MouseBox (or equivalent) configured with 2 ports.
% > Connect the left port in the box to Bpod Port#1.
%Connect the Middle Port in the box to Bpod Port #2.
% > % > Connect the right port in the box to Bpod Port#3.
% > Make sure the liquid calibration tables for ports 1 and 2 and 3 are there


global BpodSystem
global CurPokeNum
global CurRewardNum
global TrialInitialCondition
global WrongPortCount
global WrongPortThres

BpodSystem.SoftCodeHandlerFunction = 'PRSoftcodeHandler';

%% Define parameters
S = BpodSystem.ProtocolSettings; % Load settings chosen in launch manager into current workspace as a struct called S
if isempty(fieldnames(S))  % If settings file was an empty struct, populate struct with default settings
    
    S.GUI.BiasedTaskNum = 0;
    S.GUI.RewardAmountL = 2; %ul
    S.GUI.RewardAmountR = 14; %ul
    S.GUI.iti = 5; %s
    S.GUI.WaitForMidPoke = 5;
    S.GUI.TimeOut = 40; %s
    S.GUI.WrongPortThreshold = 1;
     
end

CurPokeNum = 0;
CurRewardNum = 0;
WrongPortCount = 0;

TrialInitialCondition = 0;
% 0: New session started or breakingpoint reached. Reset CurPokeNum and CurRewardNum to 0; 
% 1: Larger reward obtained. Reset only CurPokeNum to 0; 
% 2: More pokes required. Neither is reset.



% Initialize parameter GUI plugin
BpodParameterGUI('init', S);

%% Define trials
MaxTrials = 1000;
tpt       = 20;
TrialTypes = ceil(rand(1,1000)*1);
BpodSystem.Data.TrialTypes = []; % The trial type of each trial completed will be added here.

%% Initialize plots
BpodSystem.ProtocolFigures.SideOutcomePlotFig = figure('Position', [200 200 1000 200],'name','Outcome plot','numbertitle','off', 'MenuBar', 'none', 'Resize', 'off');
BpodSystem.GUIHandles.SideOutcomePlot = axes('Position', [.075 .3 .89 .6]);
SideOutcomePlot(BpodSystem.GUIHandles.SideOutcomePlot,'init',2-TrialTypes);
% BpodNotebook('init');

%% Main trial loop
if S.GUI.RewardAmountL >= S.GUI.RewardAmountR
    LargerPortIn = 'Port1In';
    SmallerPortIn = 'Port3In';
    LargerRewardValve = 1;
    SmallerRewardValve = 4;
    LargerReward = 'LeftReward';
    SmallerReward = 'RightReward';
else
    LargerPortIn = 'Port3In';
    SmallerPortIn = 'Port1In'; 
    LargerRewardValve = 4;
    SmallerRewardValve = 1;
    LargerReward = 'RightReward';
    SmallerReward = 'LeftReward';
end

for currentTrial = 1:MaxTrials
    S = BpodParameterGUI('sync', S); % Sync parameters with BpodParameterGUI plugin
    R = GetValveTimes(S.GUI.RewardAmountL, [1]);
    LeftValveTime = R(1);
    R = GetValveTimes(S.GUI.RewardAmountR, [3]);
    RightValveTime = R(1); % Update reward amounts    
    WrongPortThres = S.GUI.WrongPortThreshold;
    
    if S.GUI.RewardAmountL >= S.GUI.RewardAmountR
        R = GetValveTimes(S.GUI.RewardAmountL, [1]);
        LargerValveTime = R(1);
        R = GetValveTimes(S.GUI.RewardAmountR, [3]);
        SmallerValveTime = R(1);
    else
        R = GetValveTimes(S.GUI.RewardAmountR, [3]);
        LargerValveTime = R(1);
        R = GetValveTimes(S.GUI.RewardAmountL, [1]);
        SmallerValveTime = R(1);
    end
    
    sma = NewStateMatrix(); % Assemble state matrix
    
    if currentTrial <=  S.GUI.BiasedTaskNum
        
        sma = AddState(sma, 'Name', 'iti', ...
            'Timer', S.GUI.iti, ...
            'StateChangeConditions', {'Tup', 'WaitForMidPoke'}, ...
            'OutputActions', {});       

        sma = AddState(sma, 'Name', 'WaitForMidPoke', ...
            'Timer', tpt,...
            'StateChangeConditions', {'Tup', 'exit', 'Port2In', 'WaitForSidePoke'},...
            'OutputActions', {'PWM2', 255});

        sma = AddState(sma, 'Name', 'WaitForSidePoke', ...
            'Timer', tpt,...
            'StateChangeConditions', {'Tup', 'exit', 'Port1In', 'LeftReward', 'Port3In', 'RightReward'},...
            'OutputActions', {'PWM1', 255, 'PWM3', 255});

        sma = AddState(sma, 'Name', 'LeftReward', ...
            'Timer', LeftValveTime,...
            'StateChangeConditions', {'Tup', 'exit'},...
            'OutputActions', {'ValveState', 1});

        sma = AddState(sma, 'Name', 'RightReward', ...
            'Timer', RightValveTime,...
            'StateChangeConditions', {'Tup', 'exit'},...
            'OutputActions', {'ValveState', 4});   

    else
        
        if TrialInitialCondition == 0
        
            sma = AddState(sma, 'Name', 'iti', ...
                'Timer', S.GUI.iti, ...
                'StateChangeConditions', {'Tup', 'WaitForMidPoke'}, ...
                'OutputActions', {'SoftCode', 100});
            
            sma = AddState(sma, 'Name', 'WaitForMidPoke', ...
                'Timer', tpt,...
                'StateChangeConditions', {'Tup', 'InactiveTrial', 'Port2In', 'WaitForSidePoke'},...
                'OutputActions', {'PWM2', 255});
        
        elseif TrialInitialCondition == 1
            
            sma = AddState(sma, 'Name', 'iti', ...
                'Timer', S.GUI.iti, ...
                'StateChangeConditions', {'Tup', 'WaitForMidPoke'}, ...
                'OutputActions', {'SoftCode', 101});
            
            sma = AddState(sma, 'Name', 'WaitForMidPoke', ...
                'Timer', tpt,...
                'StateChangeConditions', {'Tup', 'InactiveTrial', 'Port2In', 'WaitForSidePoke'},...
                'OutputActions', {'PWM2', 255});

        end
        
        sma = AddState(sma, 'Name', 'WaitForSidePoke', ...
            'Timer', tpt,...
            'StateChangeConditions', {'Tup', 'InactiveTrial', LargerPortIn, 'LargerReward_Rqst', SmallerPortIn, SmallerReward},...
            'OutputActions', {'PWM1', 255, 'PWM3', 255});
        
        sma = AddState(sma, 'Name', 'LargerReward_Rqst', ...
            'Timer', 0,...
            'StateChangeConditions', {'Tup', 'LargerReward_Judge'},...
            'OutputActions', {'SoftCode', 200});

        sma = AddState(sma, 'Name', 'LargerReward_Judge', ...
            'Timer', 10,...
            'StateChangeConditions', {'Tup', 'exit', 'SoftCode1', 'exit', 'SoftCode2', LargerReward},...
            'OutputActions', {});      
        
        sma = AddState(sma, 'Name', LargerReward, ...
            'Timer', LargerValveTime,...
            'StateChangeConditions', {'Tup', 'exit'},...
            'OutputActions', {'ValveState', LargerRewardValve})  ;      
        
        sma = AddState(sma, 'Name', SmallerReward, ...
            'Timer', SmallerValveTime,...
            'StateChangeConditions', {'Tup', 'TimeOut'},...
            'OutputActions', {'ValveState', SmallerRewardValve, 'SoftCode', 150}); 
        
        sma = AddState(sma, 'Name', 'InactiveTrial',...
            'Timer', 0,...
            'StateChangeConditions', {'Tup', 'TimeOut'},...
            'OutputActions', {'SoftCode', 170});

        sma = AddState(sma, 'Name', 'TimeOut',...
            'Timer', S.GUI.TimeOut,...
            'StateChangeConditions', {'Tup', 'exit'},...
            'OutputActions', {});  
        
    end
    
    SendStateMatrix(sma);
    RawEvents = RunStateMatrix;

    disp(['CurrentTrialNum = ', num2str(currentTrial)])
    disp(['TrialInitialCondition = ' num2str(TrialInitialCondition)])
    disp(['CurPokeNum = ' num2str(CurPokeNum)])
    disp(['ProgPokeReq = ' num2str(round(5*exp(0.2*(CurRewardNum + 1))-5))])
    disp(['CurRewardNum = ' num2str(CurRewardNum)])
    disp(['WrongPortCount = ' num2str(WrongPortCount)])
    disp(['WrongPortThres = ' num2str(WrongPortThres)])
    disp('---------------')
    
    if ~isempty(fieldnames(RawEvents)) % If trial data was returned
        BpodSystem.Data = AddTrialEvents(BpodSystem.Data,RawEvents); % Computes trial events from raw data
        %BpodSystem.Data = BpodNotebook('sync', BpodSystem.Data); % Sync with Bpod notebook plugin
        BpodSystem.Data.TrialSettings(currentTrial) = S; % Adds the settings used for the current trial to the Data struct (to be saved after the trial ends)
        BpodSystem.Data.TrialInitialCondition(currentTrial) = TrialInitialCondition;
        BpodSystem.Data.CurrentPokeNumber(currentTrial) = CurPokeNum;
        BpodSystem.Data.CurrentRewardNumber(currentTrial) = CurRewardNum;
        BpodSystem.Data.ProgressivePokeRequirement(currentTrial) = round(5*exp(0.2*(CurRewardNum + 1))-5);
        BpodSystem.Data.WrongPortCount(currentTrial) = WrongPortCount;
        BpodSystem.Data.WrongPortThreshold(currentTrial) = WrongPortThres;
        UpdateSideOutcomePlot(TrialTypes, BpodSystem.Data);
        SaveBpodSessionData; % Saves the field BpodSystem.Data to the current data file
    end 
    
    HandlePauseCondition; % Checks to see if the protocol is paused. If so, waits until user resumes.
    if BpodSystem.BeingUsed == 0
        return
    end
    
end



function UpdateSideOutcomePlot(TrialTypes, Data)
   
    global BpodSystem S
    
    Outcomes = zeros(1,Data.nTrials);
    
    for x = 1:Data.nTrials
        
        if ~isnan(Data.RawEvents.Trial{x}.States.LeftReward(1))
            Outcomes(x) = 1;
        elseif ~isnan(Data.RawEvents.Trial{x}.States.RightReward(1))
            Outcomes(x) = 0;
        else
            Outcomes(x) = 3;
        end
        
    end
    
    SideOutcomePlot(BpodSystem.GUIHandles.SideOutcomePlot,'update',Data.nTrials+1,3-TrialTypes,Outcomes);

    for x = 1: length(Data.RawEvents.Trial)

        countsright(x) = ~isnan(Data.RawEvents.Trial{x}.States.RightReward(1));
        countsleft(x)  = ~isnan(Data.RawEvents.Trial{x}.States.LeftReward(1));

        if isnan(Data.RawEvents.Trial{x}.States.LeftReward(1))&&isnan(Data.RawEvents.Trial{x}.States.RightReward(1))
            timeup(x)=1;
        else
            timeup(x)=0;
        end

    end

    figure(5)

    subplot(311)
    plot(cumsum(countsright),'b'),hold on
    plot(cumsum(countsleft),'r'),hold on
    xlabel('trial number')
    ylabel('cumulative enters')
    legend('right','left')
    ylim([0 length(Data.RawEvents.Trial)])

    subplot(312)
    plot(Data.CurrentPokeNumber, 'm'),hold on
    plot(Data.ProgressivePokeRequirement, 'c'),hold on
    xlabel('trial number')
    ylabel('value')
    legend('CurPokeNum','ProgPokeReq')
%   ylim([0 length(Data.RawEvents.Trial)])

    subplot(313)
    plot(Data.CurrentRewardNumber, 'k'),hold on
    xlabel('trial number')
    ylabel('CurRewardNum')
    legend('CurRewardNum')
%   ylim([0 length(Data.RawEvents.Trial)])
    
%     subplot(414)
%     plot(Data.WrongPortCount, 'g'),hold on
%     plot(Data.WrongPortThreshold, 'y'),hold on
%     xlabel('trial number')
%     ylabel('value')
%     legend('WrongPortCount','WrongPortThres')
%     ylim([0 length(Data.RawEvents.Trial)])
