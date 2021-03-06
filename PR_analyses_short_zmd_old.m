%%% zmd mods to sarah's PR analysis code for the PR W2 cohort
%%% skipped may 31st RHRH file bc was corrupted

close all; clear all
%list of animals in cage
animals{1}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/prw2_lh';
animals{2}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/prw2_lhrh';
animals{3}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/prw2_nh';
animals{4}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/prw2_rh';
animals{5}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/prw2_rhrh';
annames = ["LH", "LHRH", "NH", "RH", "RHRH"]; %mapping functions to animal names
%dest for figures
<<<<<<< HEAD:PR_analyses_short_zmd_old.m
dst='/home/kepecs/Desktop/analysis';
=======
dst='C:\Users\zahhr\Box\kepecs_lab_summer2021\behavior_data_balbc_progressive_ratio_w2\analysis';
>>>>>>> 05432985256f9166dbad45ef31a17bb3436f8f6b:PR_analyses_short_zmd.m
back=20; %what is this for?
rel=NaN(2,5,50);
poke=NaN(2,5,50);
rewa=NaN(2,5,50);
for ani=1:5
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
        trials = {1:50, length(SessionData.RawEvents.Trial)-50+1:length(SessionData.RawEvents.Trial)}; %trail types; %only look at first 50 trials!!!!!!!!!!!
        %trials = length(SessionData.RawEvents.Trial)-30+1:length(SessionData.RawEvents.Trial);
        %length(SessionData.RawEvents.Trial) ==> for all trails
        bp=[];
        rew=[];
        for t=1:length(trials)
            trial=cell2mat(trials(t));
            for i=trial %calc breaking point
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
                end
            end
            %2 loops doing the same thing for diff vars?
            for i=trial %save reward
                if SessionData.TrialSettings(1).GUI.RewardAmountL  > SessionData.TrialSettings(1).GUI.RewardAmountR
                    if ~isnan(SessionData.RawEvents.Trial{1,i}.States.RightReward)
                        rew(t,i)=0;
                    elseif ~isnan(SessionData.RawEvents.Trial{1,i}.States.LeftReward)
                        rew(t,i)=1;
                    else
                        rew(t,i)=NaN;
                    end
                else
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
            mbp(t,f)=nanmean(bp(t,trial)); %mean breaking point for that specific session for 1 animal
            rewall{t,ani,f}=rew(t,trial); %all rewards

            si(1)=nansum(rewall{t,ani,f}==1);
            si(2)=nansum(rewall{t,ani,f}==0);
            rel(t,ani, f)=si(1)/si(2);
            rewa(t,ani,f)=si(2)*2+si(1)*14; %relative rewards
            poke(t,ani,f)=nansum(si); %number of rewards/pokes
        end
    end
    %mark day of tumor cell injection
    day = 22;
    if contains(flist(1).name, "RHRH") %because of the behav files is corrupted
        day = 21;
    end
    %format dates
    dates = {flist.date};
    for i=1:length(dates)
        tmp = char(dates(i));
        dates(i) = {tmp(1:6)};
    end

    %mark day of tumor cell injection
    day = 22;
    if contains(flist(1).name, "RHRH") %because of the behav files is corrupted
        day = 21;
    end
    %format dates
    dates = {flist.date};
    for i=1:length(dates)
        tmp = char(dates(i));
        dates(i) = {tmp(1:6)};
    end
    fig = figure();
    subplot(221);
<<<<<<< HEAD:PR_analyses_short_zmd_old.m
    plot(mbp, 'k', 'LineWidth',2), hold on
=======
    plot(mbp(1,:), 'k', 'LineWidth',2), hold on
    plot(mbp(2,:), 'm', 'LineWidth',2), hold on
>>>>>>> 05432985256f9166dbad45ef31a17bb3436f8f6b:PR_analyses_short_zmd.m
    xticks(1:length(mbp))
    xticklabels(dates)
    xline(day,'--r',{'C26','injection'});
    ylabel('mean breaking point')
    xlabel('training sessions')
<<<<<<< HEAD:PR_analyses_short_zmd_old.m
    hold off
    subplot(222)
    plot( rel(ani,:), 'k', 'LineWidth',2), hold on
=======
    legend('first 50 trials','last 50 trials')
    hold off
    subplot(222)
    plot(reshape(rel(1, ani,:), 1, []), 'k', 'LineWidth',2), hold on
    plot(reshape(rel(2, ani,:), 1, []), 'm', 'LineWidth',2), hold on
>>>>>>> 05432985256f9166dbad45ef31a17bb3436f8f6b:PR_analyses_short_zmd.m
    xline(day,'--r',{'C26','injection'});
    ylabel('relation high/low')
    xlabel('training sessions')
    xticks(1:length(mbp))
    xticklabels(dates)
<<<<<<< HEAD:PR_analyses_short_zmd_old.m
    subplot(223)
    plot( rewa(ani,:), 'k','LineWidth',2), hold on
=======
    legend('first 50 trials','last 50 trials')
    hold off
    subplot(223)
    plot(reshape(rewa(1, ani,:), 1, []), 'k', 'LineWidth',2), hold on
    plot(reshape(rewa(2, ani,:), 1, []), 'm', 'LineWidth',2), hold on
>>>>>>> 05432985256f9166dbad45ef31a17bb3436f8f6b:PR_analyses_short_zmd.m
    xline(day,'--r',{'C26','injection'});
    ylabel('water intake')
    xlabel('training sessions')
    xticks(1:length(mbp))
    xticklabels(dates)
<<<<<<< HEAD:PR_analyses_short_zmd_old.m
    subplot(224)
    plot( poke(ani,:), 'k','LineWidth',2), hold on
=======
    legend('first 50 trials','last 50 trials')
    hold off
    subplot(224)
    plot(reshape(poke(1, ani,:), 1, []), 'k', 'LineWidth',2), hold on
    plot(reshape(poke(2, ani,:), 1, []), 'm', 'LineWidth',2), hold on
>>>>>>> 05432985256f9166dbad45ef31a17bb3436f8f6b:PR_analyses_short_zmd.m
    xline(day,'--r',{'C26','injection'});
    ylabel('pokes')
    xlabel('training sessions')
    xticks(1:length(mbp))
    xticklabels(dates)
<<<<<<< HEAD:PR_analyses_short_zmd_old.m
    sgtitle(sprintf('summary for %s', annames(ani)))
    currfile = strcat(dst, '/', annames(ani), sprintf('_training_summary_sessions1-%d.jpeg', length(flist)));
    saveas(fig, currfile)
    %save to larger array for plots across animals?
    pokeall(ani,:)=poke(ani,:);
    rewaall(ani,:)=rewa(ani,:);
    relall(ani,:)=rel(ani,:);
    mbpall(ani,1:length(mbp))=mbp;
=======
    legend('first 50 trials','last 50 trials')
    sgtitle(sprintf('summary for %s', annames(ani)))
    currfile = strcat(dst, '\', annames(ani), ...
        sprintf('_training_summary_sessions1-%d_firstlast30trails.jpeg', length(flist)));
    saveas(fig, currfile)
    %save to larger array for plots across animals?
    pokeall(1, ani,:) = reshape(poke(1, ani,:), 1, []); %first 50 trials
    pokeall(2, ani,:) = reshape(poke(2, ani,:), 1, []); %last 50 trials
    rewaall(1, ani,:) = reshape(rewa(1, ani,:), 1, []);
    rewaall(2, ani,:) = reshape(rewa(2, ani,:), 1, []);
    relaall(1, ani, :) = reshape(rel(1, ani,:), 1, []);
    relaall(2, ani, :) = reshape(rel(2, ani,:), 1, []);
    mbpall(1, ani, 1:length(mbp(1, :))) = mbp(1, :);
    mbpall(2, ani, 1:length(mbp(1, :))) = mbp(2, :);
>>>>>>> 05432985256f9166dbad45ef31a17bb3436f8f6b:PR_analyses_short_zmd.m
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
    
fig = figure();
subplot(221);
plot(reshape(nanmean(mbpall(1, 1:3, :)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(nanmean(mbpall(2, 1:3, :)), 1, []), 'm', 'LineWidth',2), hold on
plot(reshape(nanmean(mbpall(1, 3:5, :)), 1, []), ':b', 'LineWidth',2), hold on
plot(reshape(nanmean(mbpall(2, 3:5, :)), 1, []), ':g', 'LineWidth',2), hold on
xticks(1:length(reshape(nanmean(mbpall(2, 3:5, :)), 1, [])))
xticklabels(dates)
xline(day,'--r',{'C26','injection'});
ylabel('mean breaking point')
xlabel('training sessions')
title('mean breaking point')
legend('tumor first 50 trials','tumor last 50 trials', 'control first 50 trials', ...
    'control last 50 trials', 'Location','northwest','NumColumns',2)
hold off
   
subplot(222)
plot(reshape(nanmean(relaall(1, 1:3, :)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(nanmean(relaall(2, 1:3, :)), 1, []), 'm', 'LineWidth',2), hold on
plot(reshape(nanmean(relaall(1, 3:5, :)), 1, []), ':b', 'LineWidth',2), hold on
plot(reshape(nanmean(relaall(2, 3:5, :)), 1, []), ':g', 'LineWidth',2), hold on
xticks(1:length(reshape(nanmean(relaall(2, 3:5, :)), 1, [])))
xticklabels(dates)
xline(day,'--r',{'C26','injection'});
ylabel('relation high/low')
<<<<<<< HEAD:PR_analyses_short_zmd_old.m
title('relation high/low rewards') 
=======
title('relation high/low rewards')
xlabel('training sessions')
legend('tumor first 50 trials','tumor last 50 trials', 'control first 50 trials', ...
    'control last 50 trials', 'Location','northwest','NumColumns',2)
hold off
>>>>>>> 05432985256f9166dbad45ef31a17bb3436f8f6b:PR_analyses_short_zmd.m

subplot(223)
plot(reshape(nanmean(rewaall(1, 1:3, :)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(nanmean(rewaall(2, 1:3, :)), 1, []), 'm', 'LineWidth',2), hold on
plot(reshape(nanmean(rewaall(1, 3:5, :)), 1, []), ':b', 'LineWidth',2), hold on
plot(reshape(nanmean(rewaall(2, 3:5, :)), 1, []), ':g', 'LineWidth',2), hold on
xticks(1:length(reshape(nanmean(rewaall(2, 3:5, :)), 1, [])))
xticklabels(dates)
xline(day,'--r',{'C26','injection'});
ylabel('water intake')
title('water intake')
xlabel('training sessions')
legend('tumor first 50 trials','tumor last 50 trials', 'control first 50 trials', ...
    'control last 50 trials', 'Location','northwest','NumColumns',2)
hold off

subplot(224)
plot(reshape(nanmean(pokeall(1, 1:3, :)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(nanmean(pokeall(2, 1:3, :)), 1, []), 'm', 'LineWidth',2), hold on
plot(reshape(nanmean(pokeall(1, 3:5, :)), 1, []), ':b', 'LineWidth',2), hold on
plot(reshape(nanmean(pokeall(2, 3:5, :)), 1, []), ':g', 'LineWidth',2), hold on
xticks(1:length(reshape(nanmean(pokeall(2, 3:5, :)), 1, [])))
xticklabels(dates)
xline(day,'--r',{'C26','injection'});
ylabel('pokes')
title('pokes')
legend('tumor first 50 trials','tumor last 50 trials', 'control first 50 trials', ...
    'control last 50 trials', 'Location','northwest','NumColumns',2)
hold off

for i=1:4
    subplot(2,2,i)
xlabel('days to endpoint')
xticks([1 14 24 34 44])
%xlim([0 17])
xticklabels({'-back', '-30', '-20','-10','0'})
end


