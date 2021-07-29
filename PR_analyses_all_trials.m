%%% zmd mods to sarah's PR analysis code for the PR cohort (s)
%%% skipped may 31st RHRH file bc was corrupted

close all; clear all
%list of animals in cage
animals{1}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/prw2_lh';
animals{2}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/prw2_lhrh';
animals{3}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/prw2_nh';
animals{4}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/prw2_rh';
animals{5}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/prw2_rhrh';
animals{6}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/PR_W1_LH';
animals{7}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/PR_W1_LHRH';
animals{8}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/PR_W1_NH';
animals{9}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/PR_W1_RH';
animals{10}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/PR_W1_RHRH';

annames = ["LHW2", "LHRHW2", "NHW2", "RHW2", "RHRHW2", ...
          "LHW1", "LHRHW1", "NHW1", "RHW1", "RHRHW1"]; %mapping functions to animal names
%dest for figures
dst='/home/kepecs/Desktop/analysis';
for ani=1:10
     cd(animals{ani}) %navigate to animal dir
     %%%sarah's way of listing sessions
     %list('fmask','*Mar*','fn','all.txt','mode','w');
     %list('fmask','*Apr*','fn','all.txt','mode','a');           
     %list('fmask','*May*','fn','all.txt','mode','a');
  
    %flist = textread(['all.txt'],'%s');
         %length(flist)
    %%%
    %find all session files in directory
    flist = dir(fullfile(animals{ani}, '*mat'));
    % SORT BY DATE
    [~,ind] = sort([flist.datenum]);
    flist = flist(ind);
    %pad mean breaking point with nans just incase session date is corrupted??
    mbp = NaN(2,5); 
    for f =1:length(flist)     
        disp(f);
        load(fullfile(animals{ani}, flist(f).name)) %load session file      
        trials = {1:length(SessionData.RawEvents.Trial)};%{length(SessionData.RawEvents.Trial)-50:length(SessionData.RawEvents.Trial)}; 
        bp=[];rew=[];rewnum=[];rxinit=[];rxside=[];inactivetime=[];
        for t=1:length(trials)
            trial=cell2mat(trials(t));
            rewnum = SessionData.CurrentRewardNumber(trial);
            for i=trial %calc breaking point
                %get reaction times
                if isfield(SessionData.RawEvents.Trial{1,i}.States, 'WaitForMidPoke')==1 && ...
                    (isfield(SessionData.RawEvents.Trial{1,i}.States, 'InactiveTrial')==0 || ...
                    isnan(SessionData.RawEvents.Trial{1,i}.States.InactiveTrial(1))) %do not want to get inactive trial rx time
                    %get difference between start and end time
                    rxinit(i) = SessionData.RawEvents.Trial{1,i}.States.WaitForMidPoke(2)-SessionData.RawEvents.Trial{1,i}.States.WaitForMidPoke(1);                    
                        
                else
                    rxinit(i) = NaN;
                end 
                if isfield(SessionData.RawEvents.Trial{1,i}.States, 'WaitForSidePoke')==1 %get difference between start and end time
                    rxside(i) = SessionData.RawEvents.Trial{1,i}.States.WaitForSidePoke(2)-SessionData.RawEvents.Trial{1,i}.States.WaitForSidePoke(1);
                else
                    rxside(i) = NaN;
                end
                %get #of inactive trials & duration of inactive trials
                if isfield(SessionData.RawEvents.Trial{1,i}.States, 'InactiveTrial')==1
                    if ~isnan(SessionData.RawEvents.Trial{1,i}.States.InactiveTrial)
                        inactivetime(i) = SessionData.RawEvents.Trial{1,i}.States.InactiveTrial(2)-SessionData.RawEvents.Trial{1,i}.States.InactiveTrial(1);
                    else
                        inactivetime(i) = NaN;
                    end
                else
                    inactivetime(i) = NaN;
                end
                if SessionData.TrialSettings(1).GUI.RewardAmountL  > SessionData.TrialSettings(1).GUI.RewardAmountR
                    if ~isnan(SessionData.RawEvents.Trial{1,i}.States.RightReward)                
                        if i==1
                            bp(t,i) = SessionData.ProgressivePokeRequirement(i);
                        else
                            bp(t,i) = SessionData.ProgressivePokeRequirement(i-1);
                        end
                    else
                        bp(t,i) =NaN;
                    end
                    if ~isnan(SessionData.RawEvents.Trial{1,i}.States.RightReward)
                        rew(t,i)=0; %save reward
                    elseif ~isnan(SessionData.RawEvents.Trial{1,i}.States.LeftReward)
                        rew(t,i)=1;
                    else
                        rew(t,i)=NaN;
                    end
                else
                    if ~isnan(SessionData.RawEvents.Trial{1,i}.States.LeftReward)
                        if i==1
                            bp(t,i) =  SessionData.ProgressivePokeRequirement(i);
                        else
                            bp(t,i) =  SessionData.ProgressivePokeRequirement(i-1)  ;
                        end
                    else
                        bp(t,i) =NaN;
                    end
                    if ~isnan(SessionData.RawEvents.Trial{1,i}.States.RightReward)
                        rew(t,i)=1;
                    elseif ~isnan(SessionData.RawEvents.Trial{1,i}.States.LeftReward)
                        rew(t,i)=0;
                    else
                        rew(t,i)=NaN;
                    end
                end
            end
            %bp(t,isnan(bp))=[]; %don't need to do this as taking nanmean
            %anyways
            rewnums{t,ani,f} = rewnum;
            rxinits(t,ani,f) = nanmean(rxinit);
            rxsides(t,ani,f) = nanmean(rxside);
            inactivesessions(t,ani,f) = sum(~isnan(inactivetime));
            mbp(t,f)=nanmean(bp(t,trial)); %mean breaking point for that specific session for 1 animal
            rewall{t,ani,f}=rew(t,trial); %all rewards

            si(1)=nansum(rewall{t,ani,f}==1);
            si(2)=nansum(rewall{t,ani,f}==0);
            rel(t,ani,f)=si(1)/si(2);
            rewa(t,ani,f)=si(2)*2+si(1)*14; %relative rewards
            poke(t,ani,f)=nansum(si); %number of rewards/pokes
            
        end
    end

    %format dates
    dates = {flist.date};
    for i=1:length(dates)
        tmp = char(dates(i));
        dates(i) = {tmp(1:6)};
    end
%     fig = figure();
%     subplot(221);
%     plot(mbp(1,:), 'k', 'LineWidth',2), hold on
%     xticks(1:length(mbp))
%     xticklabels(dates)
%     ylabel('mean breaking point')
%     xlabel('training sessions')
%     hold off
%     subplot(222)
%     plot(reshape(rel(1, ani,:), 1, []), 'k', 'LineWidth',2), hold on
%     ylabel('relation high/low')
%     xlabel('training sessions')
%     xticks(1:length(mbp))
%     xticklabels(dates)
%     hold off
%     subplot(223)
%     plot(reshape(rewa(1, ani,:), 1, []), 'k', 'LineWidth',2), hold on
%     ylabel('water intake')
%     xlabel('training sessions')
%     xticks(1:length(mbp))
%     xticklabels(dates)
%     hold off
%     subplot(224)
%     plot(reshape(poke(1, ani,:), 1, []), 'k', 'LineWidth',2), hold on
%     ylabel('pokes')
%     xlabel('training sessions')
%     xticks(1:length(mbp))
%     xticklabels(dates)
%     sgtitle(sprintf('summary for %s', annames(ani)))
%     currfile = strcat(dst, '\', annames(ani), ...
%         sprintf('_training_summary_sessions1-%d_alltrials.jpeg', length(flist)));
%     saveas(fig, currfile)
    %save to larger array for plots across animals?
    pokeall(1, ani, 1:length(poke(1, ani,:))) = reshape(poke(1, ani,:), 1, []);
    rewaall(1, ani, 1:length(rewa(1, ani,:))) = reshape(rewa(1, ani,:), 1, []);
    relaall(1, ani, 1:length(rel(1, ani,:))) = reshape(rel(1, ani,:), 1, []);
    mbpall(1, ani, 1:length(mbp(1, :))) = mbp(1, :);
end

%%
flist = dir(fullfile(animals{1}, '*mat'));
% SORT BY DATE
[~,ind] = sort([flist.datenum]);
flist = flist(ind);
%get all dates
dates = {flist.date};
for i=1:length(dates)
    tmp = char(dates(i));
    dates(i) = {tmp(1:6)};
end
%specify endpoint for animals that have reached it
%in order of animal name
%annames = ["LHW2", "LHRHW2", "NHW2", "RHW2", "RHRHW2", ...
%          "LHW1", "LHRHW1", "NHW1", "RHW1", "RHRHW1"]; %mapping functions to animal names
%may 13 RH behavior file deleted (incomplete data)
endpoint = [length(dates), length(dates), 38, length(dates), 37,...
            33,45,32,46,52];%endpoints of old cohort same as total trials...
%get measures only until endpoint and nan the rest of the days
for i=1:length(endpoint)
    mbpall(1,i,endpoint(i)+1:length(mbpall(1,i,:))) = NaN;
    relaall(1,i,endpoint(i)+1:length(relaall(1,i,:))) = NaN;
    rewaall(1,i,endpoint(i)+1:length(rewaall(1,i,:))) = NaN;
    pokeall(1,i,endpoint(i)+1:length(pokeall(1,i,:))) = NaN;
    rxsides(1,i,endpoint(i)+1:length(rxsides(1,i,:))) = NaN;
    rxinits(1,i,endpoint(i)+1:length(rxinits(1,i,:))) = NaN;
    inactivesessions(1,i,endpoint(i)+1:length(inactivesessions(1,i,:))) = NaN;
end

fig = figure();
subplot(241);
% alignt one with the lowest endpoint
%mean across groups
tumor = [reshape(mbpall(1, 1, endpoint(1)-min(endpoint)+1:endpoint(1)), 1, []);...
    reshape(mbpall(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, []);reshape(mbpall(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []);...
    reshape(mbpall(1, 8, 1:endpoint(8)), 1, []);reshape(mbpall(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, [])];
tumormean = mean(tumor,1);
tumorstd = std(tumor, 1)/sqrt(length(tumor(:,1)));
ctrl = [reshape(mbpall(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []);...
    reshape(mbpall(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, [])];
ctrlmean = mean(ctrl,1);
ctrlstd = std(ctrl,1)/sqrt(length(ctrl(:,1)));
errorbar(1:length(tumormean),tumormean, tumorstd, 'g', 'LineWidth',2), hold on
errorbar(1:length(ctrlmean),ctrlmean, ctrlstd, 'k', 'LineWidth',2), hold on
xticks([1:8:min(endpoint)-3,30,31,32])
xticklabels(["back",-23,-15,-7,-2,-1,0])
ylim([0,15])
ylabel('mean breaking point')
xlabel('days before endpoint')
title('mean breaking point')
legend('tumor', 'ctrl','Location','northwest','NumColumns',2)
hold off
   
subplot(242)
% alignt one with the lowest endpoint
%mean across groups
tumor = [reshape(relaall(1, 1, endpoint(1)-min(endpoint)+1:endpoint(1)), 1, []);...
    reshape(relaall(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, []);reshape(relaall(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []);...
    reshape(relaall(1, 8, 1:endpoint(8)), 1, []);reshape(relaall(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, [])];
tumormean = mean(tumor,1);
tumorstd = std(tumor, 1)/sqrt(length(tumor(:,1)));
ctrl = [reshape(relaall(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []);...
    reshape(relaall(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, [])];
ctrlmean = mean(ctrl,1);
ctrlstd = std(ctrl,1)/sqrt(length(ctrl(:,1)));
errorbar(1:length(tumormean),tumormean, tumorstd, 'g', 'LineWidth',2), hold on
errorbar(1:length(ctrlmean),ctrlmean, ctrlstd, 'k', 'LineWidth',2), hold on
xticks([1:8:min(endpoint)-3,30,31,32])
xticklabels(["back",-23,-15,-7,-2,-1,0])
ylabel('relation high/low')
xlabel('days before endpoint')
title('relation high/low rewards')
legend('tumor', 'ctrl','Location','northwest','NumColumns',2)
hold off

subplot(243)
% alignt one with the lowest endpoint
%mean across groups
tumor = [reshape(rewaall(1, 1, endpoint(1)-min(endpoint)+1:endpoint(1)), 1, []);...
    reshape(rewaall(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, []);reshape(rewaall(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []);...
    reshape(rewaall(1, 8, 1:endpoint(8)), 1, []);reshape(rewaall(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, [])];
tumormean = mean(tumor,1);
tumorstd = std(tumor, 1)/sqrt(length(tumor(:,1)));
ctrl = [reshape(rewaall(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []);...
    reshape(rewaall(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, [])];
ctrlmean = mean(ctrl,1);
ctrlstd = std(ctrl,1)/sqrt(length(ctrl(:,1)));
errorbar(1:length(tumormean),tumormean, tumorstd, 'g', 'LineWidth',2), hold on
errorbar(1:length(ctrlmean),ctrlmean, ctrlstd, 'k', 'LineWidth',2), hold on
xticks([1:8:min(endpoint)-3,30,31,32])
xticklabels(["back",-23,-15,-7,-2,-1,0])
ylabel('water intake')
title('water intake')
xlabel('days before endpoint')
legend('tumor', 'ctrl','Location','northwest','NumColumns',2)
hold off

subplot(244)
% alignt one with the lowest endpoint
%mean across groups
tumor = [reshape(pokeall(1, 1, endpoint(1)-min(endpoint)+1:endpoint(1)), 1, []);...
    reshape(pokeall(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, []);reshape(pokeall(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []);...
    reshape(pokeall(1, 8, 1:endpoint(8)), 1, []);reshape(pokeall(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, [])];
tumormean = mean(tumor,1);
tumorstd = std(tumor, 1)/sqrt(length(tumor(:,1)));
ctrl = [reshape(pokeall(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []);...
    reshape(pokeall(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, [])];
ctrlmean = mean(ctrl,1);
ctrlstd = std(ctrl,1)/sqrt(length(ctrl(:,1)));
errorbar(1:length(tumormean),tumormean, tumorstd, 'g', 'LineWidth',2), hold on
errorbar(1:length(ctrlmean),ctrlmean, ctrlstd, 'k', 'LineWidth',2), hold on
xticks([1:8:min(endpoint)-3,30,31,32])
xticklabels(["back",-23,-15,-7,-2,-1,0])
ylabel('pokes')
title('pokes')
xlabel('days before endpoint')
legend('tumor', 'ctrl','Location','northwest','NumColumns',2)
hold off

subplot(245)
% alignt one with the lowest endpoint
%mean across groups
tumor = [reshape(rxinits(1, 1, endpoint(1)-min(endpoint)+1:endpoint(1)), 1, []);...
    reshape(rxinits(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, []);reshape(rxinits(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []);...
    reshape(rxinits(1, 8, 1:endpoint(8)), 1, []);reshape(rxinits(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, [])];
tumormean = mean(tumor,1);
tumorstd = std(tumor, 1)/sqrt(length(tumor(:,1)));
ctrl = [reshape(rxinits(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []);...
    reshape(rxinits(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, [])];
ctrlmean = mean(ctrl,1);
ctrlstd = std(ctrl,1)/sqrt(length(ctrl(:,1)));
errorbar(1:length(tumormean),tumormean, tumorstd, 'g', 'LineWidth',2), hold on
errorbar(1:length(ctrlmean),ctrlmean, ctrlstd, 'k', 'LineWidth',2), hold on
xticks([1:8:min(endpoint)-3,30,31,32])
xticklabels(["back",-23,-15,-7,-2,-1,0])
ylabel('mean reaction time (s)')
xlabel('days before endpoint')
title('reaction time to initialize trial')
legend('tumor', 'ctrl','Location','northwest','NumColumns',2)
hold off

subplot(246)
% alignt one with the lowest endpoint
%mean across groups
tumor = [reshape(rxsides(1, 1, endpoint(1)-min(endpoint)+1:endpoint(1)), 1, []);...
    reshape(rxsides(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, []);reshape(rxsides(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []);...
    reshape(rxsides(1, 8, 1:endpoint(8)), 1, []);reshape(rxsides(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, [])];
tumormean = mean(tumor,1);
tumorstd = std(tumor, 1)/sqrt(length(tumor(:,1)));
ctrl = [reshape(rxsides(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []);...
    reshape(rxsides(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, [])];
ctrlmean = mean(ctrl,1);
ctrlstd = std(ctrl,1)/sqrt(length(ctrl(:,1)));
errorbar(1:length(tumormean),tumormean, tumorstd, 'g', 'LineWidth',2), hold on
errorbar(1:length(ctrlmean),ctrlmean, ctrlstd, 'k', 'LineWidth',2), hold on
xticks([1:8:min(endpoint)-3,30,31,32])
xticklabels(["back",-23,-15,-7,-2,-1,0])
ylabel('mean reaction time (s)')
title('reaction time to poke port')
xlabel('days before endpoint')
legend('tumor', 'ctrl','Location','northwest','NumColumns',2)

subplot(247)
% alignt one with the lowest endpoint
%mean across groups
tumor = [reshape(inactivesessions(1, 1, endpoint(1)-min(endpoint)+1:endpoint(1)), 1, []);...
    reshape(inactivesessions(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, []);reshape(inactivesessions(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []);...
    reshape(inactivesessions(1, 8, 1:endpoint(8)), 1, []);reshape(inactivesessions(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, [])];
tumormean = mean(tumor,1);
tumorstd = std(tumor, 1)/sqrt(length(tumor(:,1)));
ctrl = [reshape(inactivesessions(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []);...
    reshape(inactivesessions(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, [])];
ctrlmean = mean(ctrl,1);
ctrlstd = std(ctrl,1)/sqrt(length(ctrl(:,1)));
errorbar(1:length(tumormean),tumormean, tumorstd, 'g', 'LineWidth',2), hold on
errorbar(1:length(ctrlmean),ctrlmean, ctrlstd, 'k', 'LineWidth',2), hold on
xticks([1:8:min(endpoint)-3,30,31,32])
xticklabels(["back",-23,-15,-7,-2,-1,0])
ylabel('# of inactive trials')
title('total number of inactive trials per day')
xlabel('days before endpoint')
legend('tumor', 'ctrl','Location','northwest','NumColumns',2)
ax = subplot(248);
text(0,0.5,'tumor group, n = 5, control group n = 2');
text(0,0.3, 'error bars represent standard error');
set(ax,'visible','off')
sgtitle('behavior variables for progressive ratio task')
currfile = strcat(dst, '/', "behavior_variables_summary_prw1_prw2.fig");
saveas(fig, currfile)
%%
subplot(243)
% alignt one with the lowest endpoint
plot(reshape(rewaall(1, 1, endpoint(1)-min(endpoint)+1:endpoint(1)), 1, []), 'g', 'LineWidth',2), hold on
plot(reshape(rewaall(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, []), 'g', 'LineWidth',2), hold on
plot(reshape(rewaall(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(rewaall(1, 8, :), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(rewaall(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(rewaall(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(rewaall(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(rewaall(1, 7, endpoint(7)-min(endpoint)+1:endpoint(7)), 1, []), '--k', 'LineWidth',2), hold on
plot(reshape(rewaall(1, 9, endpoint(9)-min(endpoint)+1:endpoint(9)), 1, []), '--k', 'LineWidth',2), hold on
xticks([1:8:min(endpoint),32])
xticklabels(["back",-24,-16,-8,0])
ylabel('water intake')
title('water intake')
xlabel('days before endpoint')
legend('LHW2', 'NHW2', 'LHW1', 'NHW1', 'RHRHW1', ...
    'RHRHW2', 'RHW2', 'LHRHW1', 'RHW1', 'Location','northwest','NumColumns',2)
hold off

subplot(244)
% alignt one with the lowest endpoint
plot(reshape(pokeall(1, 1, endpoint(1)-min(endpoint)+1:endpoint(1)), 1, []), 'g', 'LineWidth',2), hold on
plot(reshape(pokeall(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, []), 'g', 'LineWidth',2), hold on
plot(reshape(pokeall(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(pokeall(1, 8, :), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(pokeall(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(pokeall(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(pokeall(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(pokeall(1, 7, endpoint(7)-min(endpoint)+1:endpoint(7)), 1, []), '--k', 'LineWidth',2), hold on
plot(reshape(pokeall(1, 9, endpoint(9)-min(endpoint)+1:endpoint(9)), 1, []), '--k', 'LineWidth',2), hold on
xticks([1:8:min(endpoint),32])
xticklabels(["back",-24,-16,-8,0])
ylabel('pokes')
title('pokes')
xlabel('days before endpoint')
legend('LHW2', 'NHW2', 'LHW1', 'NHW1', 'RHRHW1', ...
    'RHRHW2', 'RHW2', 'LHRHW1', 'RHW1', 'Location','northwest','NumColumns',2)
hold off

subplot(245)
% alignt one with the lowest endpoint
plot(reshape(rxinits(1, 1, endpoint(1)-min(endpoint)+1:endpoint(1)), 1, []), 'g', 'LineWidth',2), hold on
plot(reshape(rxinits(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, []), 'g', 'LineWidth',2), hold on
plot(reshape(rxinits(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(rxinits(1, 8, :), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(rxinits(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(rxinits(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(rxinits(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(rxinits(1, 7, endpoint(7)-min(endpoint)+1:endpoint(7)), 1, []), '--k', 'LineWidth',2), hold on
plot(reshape(rxinits(1, 9, endpoint(9)-min(endpoint)+1:endpoint(9)), 1, []), '--k', 'LineWidth',2), hold on
xticks([1:8:min(endpoint),32])
xticklabels(["back",-24,-16,-8,0])
ylabel('mean reaction time (s)')
title('reaction time to initialize trial')
xlabel('days before endpoint')
legend('LHW2', 'NHW2', 'LHW1', 'NHW1', 'RHRHW1', ...
    'RHRHW2', 'RHW2', 'LHRHW1', 'RHW1', 'Location','northwest','NumColumns',2)
hold off

subplot(246)
% alignt one with the lowest endpoint
plot(reshape(rxsides(1, 1, endpoint(1)-min(endpoint)+1:endpoint(1)), 1, []), 'g', 'LineWidth',2), hold on
plot(reshape(rxsides(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, []), 'g', 'LineWidth',2), hold on
plot(reshape(rxsides(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(rxsides(1, 8, :), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(rxsides(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(rxsides(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(rxsides(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(rxsides(1, 7, endpoint(7)-min(endpoint)+1:endpoint(7)), 1, []), '--k', 'LineWidth',2), hold on
plot(reshape(rxsides(1, 9, endpoint(9)-min(endpoint)+1:endpoint(9)), 1, []), '--k', 'LineWidth',2), hold on
xticks([1:8:min(endpoint),32])
xticklabels(["back",-24,-16,-8,0])
ylabel('mean reaction time (s)')
title('reaction time to poke port')
xlabel('days before endpoint')
legend('LHW2', 'NHW2', 'LHW1', 'NHW1', 'RHRHW1', ...
    'RHRHW2', 'RHW2', 'LHRHW1', 'RHW1', 'Location','northwest','NumColumns',2)
hold off

subplot(247)
% alignt one with the lowest endpoint
plot(reshape(inactivesessions(1, 1, endpoint(1)-min(endpoint)+1:endpoint(1)), 1, []), 'g', 'LineWidth',2), hold on
plot(reshape(inactivesessions(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, []), 'g', 'LineWidth',2), hold on
plot(reshape(inactivesessions(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(inactivesessions(1, 8, :), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(inactivesessions(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(inactivesessions(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(inactivesessions(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(inactivesessions(1, 7, endpoint(7)-min(endpoint)+1:endpoint(7)), 1, []), '--k', 'LineWidth',2), hold on
plot(reshape(inactivesessions(1, 9, endpoint(9)-min(endpoint)+1:endpoint(9)), 1, []), '--k', 'LineWidth',2), hold on
xticks([1:8:min(endpoint),32])
xticklabels(["back",-24,-16,-8,0])
ylabel('# of inactive sessions')
title('total number of inactive sessions per day')
xlabel('days before endpoint')
legend('LHW2', 'NHW2', 'LHW1', 'NHW1', 'RHRHW1', ...
    'RHRHW2', 'RHW2', 'LHRHW1', 'RHW1', 'Location','northwest','NumColumns',2)
hold off

%%
%plot trial behavior at last day (endpoint) for all animals
fig = figure();
subplot(331)
plot(cell2mat(rewnums(1,1,endpoint(1))),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,1,endpoint(1)-1)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,1,endpoint(1)-2)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,1,endpoint(1)-20)),'k','LineWidth',2),hold on
title('trial activity for LHW2')
xlabel('trial number')
ylabel('CurRewardNum')
legend('last day', 'last day-1', 'last day-2','last day-20')
subplot(332)
plot(cell2mat(rewnums(1,3,endpoint(3))),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,3,endpoint(3)-1)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,3,endpoint(3)-2)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,3,endpoint(3)-20)),'k','LineWidth',2),hold on
title('trial activity for NHW2')
xlabel('trial number')
ylabel('CurRewardNum')
legend('last day', 'last day-1', 'last day-2','last day-20')
ylabel('CurRewardNum')
subplot(333)
plot(cell2mat(rewnums(1,6,endpoint(6))),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,6,endpoint(6)-1)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,6,endpoint(6)-2)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,6,endpoint(6)-20)),'k','LineWidth',2),hold on
title('trial activity for LHW1')
xlabel('trial number')
ylabel('CurRewardNum')
legend('last day', 'last day-1', 'last day-2', 'last day-20')
subplot(334)
plot(cell2mat(rewnums(1,8,endpoint(8))),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,8,endpoint(8)-1)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,8,endpoint(8)-2)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,8,endpoint(8)-20)),'k','LineWidth',2),hold on
title('trial activity for NHW1')
xlabel('trial number')
ylabel('CurRewardNum')
legend('last day', 'last day-1', 'last day-2', 'last day-20')
subplot(335)
plot(cell2mat(rewnums(1,10,endpoint(10))),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,10,endpoint(10)-1)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,10,endpoint(10)-2)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,10,endpoint(10)-20)),'k','LineWidth',2),hold on
title('trial activity for RHRHW1')
xlabel('trial number')
ylabel('CurRewardNum')
legend('last day', 'last day-1', 'last day-2', 'last day-20')
subplot(336)
plot(cell2mat(rewnums(1,4,endpoint(4))),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,4,endpoint(4)-1)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,4,endpoint(4)-2)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,4,endpoint(4)-20)),'k','LineWidth',2),hold on
title('trial activity for RHRHW2 at last days')
xlabel('trial number')
ylabel('CurRewardNum')
legend('last day', 'last day-1', 'last day-2', 'last day-20')
subplot(337)
plot(cell2mat(rewnums(1,2,endpoint(2))),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,2,endpoint(2)-1)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,2,endpoint(2)-2)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,2,endpoint(2)-20)),'k','LineWidth',2),hold on
title('trial activity for RHW2')
xlabel('trial number')
ylabel('CurRewardNum')
legend('last day', 'last day-1', 'last day-2', 'last day-20')
subplot(338)
plot(cell2mat(rewnums(1,7,endpoint(7))),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,7,endpoint(7)-1)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,7,endpoint(7)-2)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,7,endpoint(7)-20)),'k','LineWidth',2),hold on
title('trial activity for LHRHW1')
xlabel('trial number')
ylabel('CurRewardNum')
legend('last day', 'last day-1', 'last day-2', 'last day-20')
subplot(339)
plot(cell2mat(rewnums(1,9,endpoint(9))),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,9,endpoint(9)-1)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,9,endpoint(9)-2)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,9,endpoint(9)-22)),'k','LineWidth',2),hold on
title('trial activity for RHW1')
xlabel('trial number')
ylabel('CurRewardNum')
legend('last day', 'last day-1', 'last day-2', 'last day-22')
% plot(cell2mat(rewnums(1,6,endpoint(6))),'--g','LineWidth',2),hold on
% plot(cell2mat(rewnums(1,8,endpoint(8))),'--g','LineWidth',2),hold on
% plot(cell2mat(rewnums(1,10,endpoint(10))),'--g','LineWidth',2),hold on
% plot(cell2mat(rewnums(1,7,endpoint(7))),'--k','LineWidth',2),hold on
% plot(cell2mat(rewnums(1,9,endpoint(9))),'--k','LineWidth',2),hold on

%%
%to plot per animal
plot(reshape(rewaall(1, 1, endpoint(1)-min(endpoint)+1:endpoint(1)), 1, []), 'g', 'LineWidth',2), hold on
plot(reshape(rewaall(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, []), 'g', 'LineWidth',2), hold on
plot(reshape(rewaall(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(rewaall(1, 8, :), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(rewaall(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(rewaall(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(rewaall(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(rewaall(1, 7, endpoint(7)-min(endpoint)+1:endpoint(7)), 1, []), '--k', 'LineWidth',2), hold on
plot(reshape(rewaall(1, 9, endpoint(9)-min(endpoint)+1:endpoint(9)), 1, []), '--k', 'LineWidth',2), hold on
