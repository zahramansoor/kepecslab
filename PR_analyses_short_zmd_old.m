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
dst='/home/kepecs/Desktop/analysis';
back=20; %what is this for?
rel=NaN(5,50);
poke=NaN(5,50);
rewa=NaN(5,50);
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
    mbp = NaN(1,5);
    for f =1:length(flist)     
        disp(f);
        bp=[];
        rew=[];
        load(fullfile(animals{ani}, flist(f).name)) %load session file      
        for i=1:length(SessionData.RawEvents.Trial) %calc breaking point
            if SessionData.TrialSettings(1).GUI.RewardAmountL  > SessionData.TrialSettings(1).GUI.RewardAmountR
                if ~isnan(SessionData.RawEvents.Trial{1,i}.States.RightReward)                
                    if i==1
                        bp(i) = SessionData.ProgressivePokeRequirement(i);
                    else
                        bp(i) = SessionData.ProgressivePokeRequirement(i-1);
                    end
                else
                    bp(i) =NaN;
                end
            else
                if ~isnan(SessionData.RawEvents.Trial{1,i}.States.LeftReward)
                    if i==1
                        bp(i) =  SessionData.ProgressivePokeRequirement(i);
                    else
                        bp(i) =  SessionData.ProgressivePokeRequirement(i-1)  ;
                    end
                else
                    bp(i) =NaN;
                end
            end
        end
        %2 loops doing the same thing for diff vars?
        for i=1:length(SessionData.RawEvents.Trial) %save reward
            if SessionData.TrialSettings(1).GUI.RewardAmountL  > SessionData.TrialSettings(1).GUI.RewardAmountR
                if ~isnan(SessionData.RawEvents.Trial{1,i}.States.RightReward)
                    rew(i)=0;
                elseif ~isnan(SessionData.RawEvents.Trial{1,i}.States.LeftReward)
                    rew(i)=1;
                else
                    rew(i)=NaN;
                end
            else
                if ~isnan(SessionData.RawEvents.Trial{1,i}.States.RightReward)
                    
                    rew(i)=1;
                    
                elseif ~isnan(SessionData.RawEvents.Trial{1,i}.States.LeftReward)
                    rew(i)=0;
                else
                    rew(i)=NaN;
                end
                
            end
        end
        
        bp(isnan(bp))=[]; %breaking point
        mbp(f)=mean(bp); %mean breaking point for that specific session for 1 animal
        rewall{ani,f}=rew; %all rewards
        
        si(1)=sum(rewall{ani,f}==1);
        si(2)=sum(rewall{ani,f}==0);
        rel(ani, f)=si(1)/si(2);
        rewa(ani,f)=si(2)*2+si(1)*14; %relative rewards
        poke(ani,f)=sum(si); %number of rewards/pokes
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
    plot(mbp, 'k', 'LineWidth',2), hold on
    xticks(1:length(mbp))
    xticklabels(dates)
    xline(day,'--r',{'C26','injection'});
    ylabel('mean breaking point')
    xlabel('training sessions')
    hold off
    subplot(222)
    plot( rel(ani,:), 'k', 'LineWidth',2), hold on
    xline(day,'--r',{'C26','injection'});
    ylabel('relation high/low')
    xlabel('training sessions')
    xticks(1:length(mbp))
    xticklabels(dates)
    subplot(223)
    plot( rewa(ani,:), 'k','LineWidth',2), hold on
    xline(day,'--r',{'C26','injection'});
    ylabel('water intake')
    xlabel('training sessions')
    xticks(1:length(mbp))
    xticklabels(dates)
    subplot(224)
    plot( poke(ani,:), 'k','LineWidth',2), hold on
    xline(day,'--r',{'C26','injection'});
    ylabel('pokes')
    xlabel('training sessions')
    xticks(1:length(mbp))
    xticklabels(dates)
    sgtitle(sprintf('summary for %s', annames(ani)))
    currfile = strcat(dst, '/', annames(ani), sprintf('_training_summary_sessions1-%d.jpeg', length(flist)));
    saveas(fig, currfile)
    %save to larger array for plots across animals?
    pokeall(ani,:)=poke(ani,:);
    rewaall(ani,:)=rewa(ani,:);
    relall(ani,:)=rel(ani,:);
    mbpall(ani,1:length(mbp))=mbp;
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



