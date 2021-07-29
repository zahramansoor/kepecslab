%%% zmd mods to sarah's PR analysis code for the PR W2 cohort
%%% skipped may 31st RHRH file bc was corrupted

close all; clear all
%list of animals in cage
%%% zmd mods to sarah's PR analysis code for the PR W2 cohort
%%% skipped may 31st RHRH file bc was corrupted

close all; clear all
%list of animals in cage
animals{1}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/prw2_lh';
animals{2}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/prw2_lhrh';
animals{3}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/prw2_nh';
animals{4}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/prw2_rh';
animals{5}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/prw2_rhrh';
animals{6}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/PR_W1_LH';
animals{5}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/PR_W1_LHRH';
animals{7}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/PR_W1_NH';
animals{5}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/PR_W1_RH';
animals{5}='/home/kepecs/Desktop/behavior_data_balbc_progressive_ratio_w2/PR_W1_RHRH';

annames = ["LHW2", "LHRHW2", "NHW2", "RHW2", "RHRHW2", ...
          "LHW1", "LHRHW1", "NHW1", "RHW1", "RHRHW1"]; %mapping functions to animal names
%dest for figures
dst='/home/kepecs/Desktop/analysis';
rel=NaN(2,10,50);
poke=NaN(2,10,50);
rewa=NaN(2,10,50);
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
        trials = {length(SessionData.RawEvents.Trial), length(SessionData.RawEvents.Trial)-30+1:length(SessionData.RawEvents.Trial)}; %trail types; %only look at first 50 trials!!!!!!!!!!!
        %trials = length(SessionData.RawEvents.Trial)-30+1:length(SessionData.RawEvents.Trial);
        %length(SessionData.RawEvents.Trial) ==> for all trails
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
            mbp(t,f)=nanmean(bp); %mean breaking point for that specific session for 1 animal
            rewall{t,ani,f}=bp; %all rewards

            si(1)=sum(rewall{t,ani,f}==1);
            si(2)=sum(rewall{t,ani,f}==0);
            rel(t,ani, f)=si(1)/si(2);
            rewa(t,ani,f)=si(2)*2+si(1)*14; %relative rewards
            poke(t,ani,f)=sum(si); %number of rewards/pokes
        end
    end

    %format dates
    dates = {flist.date};
    for i=1:length(dates)
        tmp = char(dates(i));
        dates(i) = {tmp(1:6)};
    end
    fig = figure();
    subplot(221);
    plot(mbp(1,:), 'k', 'LineWidth',2), hold on
    plot(mbp(2,:), 'm', LineWidth',2, 'm'), hold on
    xticks(1:length(mbp))
    xticklabels(dates)
    ylabel('mean breaking point')
    xlabel('training sessions')
    hold off
    subplot(222)
    plot( rel(1, ani,:), 'k', 'LineWidth',2), hold on
    plot( rel(2, ani,:), 'm', 'LineWidth',2), hold on
    xline(day,'--r',{'C26','injection'});
    ylabel('relation high/low')
    xlabel('training sessions')
    xticks(1:length(mbp))
    xticklabels(dates)
    subplot(223)
    plot( rewa(1,ani,:), 'k', 'LineWidth',2), hold on
    plot( rewa(2,ani,:), 'm', 'LineWidth',2), hold on
    xline(day,'--r',{'C26','injection'});
    ylabel('water intake')
    xlabel('training sessions')
    xticks(1:length(mbp))
    xticklabels(dates)
    subplot(224)
    plot( poke(1, ani,:), 'k', 'LineWidth',2), hold on
    plot( poke(2, ani,:), 'm', 'LineWidth',2), hold on
    xline(day,'--r',{'C26','injection'});
    ylabel('pokes')
    xlabel('training sessions')
    xticks(1:length(mbp))
    xticklabels(dates)
    sgtitle(sprintf('summary for %s', annames(ani)))
    currfile = strcat(dst, '\', annames(ani), ...
        sprintf('_training_summary_sessions1-%d_firstlast30trails.jpeg', length(flist)));
    saveas(fig, currfile)
    %save to larger array for plots across animals?
    %pokeall(ani,:)=poke(ani,:);
    %rewaall(ani,:)=rewa(ani,:);
    %relall(ani,:)=rel(ani,:);
    %mbpall(ani,:)=mbp;
end

%%

subplot(221)
legend('tumor','control')
ylabel('mean breaking point')
title('breaking point')
legend('LHRH','NH','LH','RHRH','RH')

subplot(222)
ylabel('relation high/low')
title('relation high/low rewards')

subplot(223)
ylabel('waterintake')
title('waterintake')

subplot(224)
ylabel('rewards received')
title('rewards received')

for i=1:4
    subplot(2,2,i)
xlabel('days to endpoint')
xticks([1 6 11 16])
%xlim([0 17])
xticklabels({'-back','-10','-5','0'})
end

%%

annames = ["LH", "LHRH", "NH", "RH", "RHRH"]; %mapping functions to animal names
%dest for figures
dst='C:\Users\zahhr\Box\kepecs_lab_summer2021\behavior_data_balbc_progressive_ratio_w2\analysis';
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
        trials = {1:30, length(SessionData.RawEvents.Trial)-30+1:length(SessionData.RawEvents.Trial)}; %trail types; %only look at first 50 trials!!!!!!!!!!!
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
            rewall{t,ani,f}=bp(t,trial); %all rewards

            si(1)=sum(rewall{t,ani,f}==1);
            si(2)=sum(rewall{t,ani,f}==0);
            rel(t,ani, f)=si(1)/si(2);
            rewa(t,ani,f)=si(2)*2+si(1)*14; %relative rewards
            poke(t,ani,f)=sum(si); %number of rewards/pokes
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
    fig = figure();
    subplot(221);
    plot(mbp(1,:), 'LineWidth',2, 'k'), hold on
    plot(mbp(2,:), 'LineWidth',2, 'm'), hold on
    xticks(1:length(mbp))
    xticklabels(dates)
    xline(day,'--r',{'C26','injection'});
    ylabel('mean breaking point')
    xlabel('training sessions')
    hold off
    subplot(222)
    plot( rel(1, ani,:), 'LineWidth',2, 'k'), hold on
    plot( rel(2, ani,:), 'LineWidth',2, 'm'), hold on
    xline(day,'--r',{'C26','injection'});
    ylabel('relation high/low')
    xlabel('training sessions')
    xticks(1:length(mbp))
    xticklabels(dates)
    subplot(223)
    plot( rewa(1,ani,:), 'LineWidth',2, 'k'), hold on
    plot( rewa(2,ani,:), 'LineWidth',2, 'm'), hold on
    xline(day,'--r',{'C26','injection'});
    ylabel('water intake')
    xlabel('training sessions')
    xticks(1:length(mbp))
    xticklabels(dates)
    subplot(224)
    plot( poke(1, ani,:), 'LineWidth',2, 'k'), hold on
    plot( poke(2, ani,:), 'LineWidth',2, 'm'), hold on
    xline(day,'--r',{'C26','injection'});
    ylabel('pokes')
    xlabel('training sessions')
    xticks(1:length(mbp))
    xticklabels(dates)
    sgtitle(sprintf('summary for %s', annames(ani)))
    currfile = strcat(dst, '\', annames(ani), ...
        sprintf('_training_summary_sessions1-%d_firstlast30trails.jpeg', length(flist)));
    saveas(fig, currfile)
    %save to larger array for plots across animals?
    %pokeall(ani,:)=poke(ani,:);
    %rewaall(ani,:)=rewa(ani,:);
    %relall(ani,:)=rel(ani,:);
    %mbpall(ani,:)=mbp;
end

%%

subplot(221)
legend('tumor','control')
ylabel('mean breaking point')
title('breaking point')
legend('LHRH','NH','LH','RHRH','RH')

subplot(222)
ylabel('relation high/low')
title('relation high/low rewards')

subplot(223)
ylabel('waterintake')
title('waterintake')

subplot(224)
ylabel('rewards received')
title('rewards received')

for i=1:4
    subplot(2,2,i)
xlabel('days to endpoint')
xticks([1 6 11 16])
%xlim([0 17])
xticklabels({'-back','-10','-5','0'})
end

%%
