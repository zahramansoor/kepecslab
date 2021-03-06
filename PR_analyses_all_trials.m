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
            %bp(isnan(bp))=[]; %don't need to do this as taking nanmean
            %anyways         
            %collect bps for cdf
            bps{t,ani,f} = bp;
            rewnums{t,ani,f} = rewnum;
            rxinits(t,ani,f) = nanmean(rxinit);
            rxsides(t,ani,f) = nanmean(rxside);
            inactivesessions(t,ani,f) = sum(~isnan(inactivetime));
            mbp(t,f)=mean(bp(t,trial)); %mean breaking point for that specific session for 1 animal
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
    bpsall(1, ani, 1:length(bps(1, ani, :))) = reshape(bps(1,ani,:),1,[]);
end

%%
%specify endpoint for animals that have reached it
%in order of animal name
%annames = ["LHW2", "LHRHW2", "NHW2", "RHW2", "RHRHW2", ...
%          "LHW1", "LHRHW1", "NHW1", "RHW1", "RHRHW1"]; %mapping functions to animal names
%old cohort = NH, LH and RHRH were cachectic
%may 13 w1RH behavior file deleted (incomplete data)
endpoint = [47, 53, 38, 53, 37,...
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
%figure for cdfs of bps for cachexia and control mice
%idea is to see whether the distribution of breaking point changes
fig1 = figure();
subplot(221);
%tumor cdfs
bpcdf = cell2mat(bpsall(1,1,endpoint(1)));
bpcdf(isnan(bpcdf))=[];
[N1,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N1),'g','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,2,endpoint(2)));
bpcdf(isnan(bpcdf))=[];
[N2,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N2),':g','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,3,endpoint(3)));
bpcdf(isnan(bpcdf))=[];
[N3,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N3),'--g','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,6,endpoint(6)));
bpcdf(isnan(bpcdf))=[];
[N6,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N6),'y','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,8,endpoint(8)));
bpcdf(isnan(bpcdf))=[];
[N8,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N8),':y','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,10,endpoint(10)));
bpcdf(isnan(bpcdf))=[];
[N10,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N10),'--y','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,4,endpoint(4)));
bpcdf(isnan(bpcdf))=[];
[N4,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N4),'k','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,5,endpoint(5)));
bpcdf(isnan(bpcdf))=[];
[N5,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N5),':k','LineWidth',2)
legend("LHW2 last day","LHRHW2 last day","NHW2 last day", "LHW1 last day", ...
    "NHW1 last day", "RHRHW1 last day","RHW2 last day", "RHRHW2 last day");
title('cdf of breaking points on last day')
xlim([0,40])
xlabel('breaking point')
ylabel('probability')

subplot(222);
bpcdf = cell2mat(bpsall(1,1,endpoint(1)-2));
bpcdf(isnan(bpcdf))=[];
[N1,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N1),'g','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,2,endpoint(2)-2));
bpcdf(isnan(bpcdf))=[];
[N2,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N2),':g','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,3,endpoint(3)-2));
bpcdf(isnan(bpcdf))=[];
[N3,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N3),'--g','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,6,endpoint(6)-2));
bpcdf(isnan(bpcdf))=[];
[N6,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N6),'y','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,8,endpoint(8)-2));
bpcdf(isnan(bpcdf))=[];
[N8,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N8),':y','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,10,endpoint(10)-2));
bpcdf(isnan(bpcdf))=[];
[N10,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N10),'--y','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,4,endpoint(4)-2));
bpcdf(isnan(bpcdf))=[];
[N4,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N4),'k','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,5,endpoint(5)-2));
bpcdf(isnan(bpcdf))=[];
[N5,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N5),':k','LineWidth',2)
legend("LHW2 last day-2","LHRHW2 last day-2","NHW2 last day-2", "LHW1 last day-2", ...
    "NHW1 last day-2", "RHRHW1 last day-2","RHW2 last day-2", "RHRHW2 last day-2");
title('cdf of breaking points on last day-2')
xlim([0,40])
xlabel('breaking point')
ylabel('probability')

subplot(223);
bpcdf = cell2mat(bpsall(1,1,endpoint(1)-15));
bpcdf(isnan(bpcdf))=[];
[N1,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N1),'g','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,2,endpoint(2)-15));
bpcdf(isnan(bpcdf))=[];
[N2,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N2),':g','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,3,endpoint(3)-15));
bpcdf(isnan(bpcdf))=[];
[N3,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N3),'--g','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,6,endpoint(6)-15));
bpcdf(isnan(bpcdf))=[];
[N6,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N6),'y','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,8,endpoint(8)-15));
bpcdf(isnan(bpcdf))=[];
[N8,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N8),':y','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,10,endpoint(10)-15));
bpcdf(isnan(bpcdf))=[];
[N10,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N10),'--y','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,4,endpoint(4)-15));
bpcdf(isnan(bpcdf))=[];
[N4,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N4),'k','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,5,endpoint(5)-15));
bpcdf(isnan(bpcdf))=[];
[N5,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N5),':k','LineWidth',2)
legend("LHW2 last day-15","LHRHW2 last day-15","NHW2 last day-15", "LHW1 last day-15", ...
    "NHW1 last day-15", "RHRHW1 last day-15","RHW2 last day-15", "RHRHW2 last day-15");
title('cdf of breaking points on last day-15')
xlim([0,40])
xlabel('breaking point')
ylabel('probability')

subplot(224);
bpcdf = cell2mat(bpsall(1,1,endpoint(1)-30));
bpcdf(isnan(bpcdf))=[];
[N1,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N1),'g','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,2,endpoint(2)-30));
bpcdf(isnan(bpcdf))=[];
[N2,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N2),':g','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,3,endpoint(3)-30));
bpcdf(isnan(bpcdf))=[];
[N3,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N3),'--g','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,6,endpoint(6)-30));
bpcdf(isnan(bpcdf))=[];
[N6,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N6),'y','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,8,endpoint(8)-30));
bpcdf(isnan(bpcdf))=[];
[N8,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N8),':y','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,10,endpoint(10)-30));
bpcdf(isnan(bpcdf))=[];
[N10,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N10),'--y','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,4,endpoint(4)-30));
bpcdf(isnan(bpcdf))=[];
[N4,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N4),'k','LineWidth',2); hold on
bpcdf = cell2mat(bpsall(1,5,endpoint(5)-30));
bpcdf(isnan(bpcdf))=[];
[N5,~] = histcounts(bpcdf,'Normalization','probability');
plot(cumsum(N5),':k','LineWidth',2)
legend("LHW2 last day-30","LHRHW2 last day-30","NHW2 last day-30", "LHW1 last day-30", ...
    "NHW1 last day-30", "RHRHW1 last day-30","RHW2 last day-30", "RHRHW2 last day-30");
title('cdf of breaking points on last day-30 (pre-tumor)')
xlim([0,40])
xlabel('breaking point')
ylabel('probability')
%%
%take mean across group, and plot with time in different colors
%endpoint
clear bpprobs;
%tumor group
tumorind = [1 2 3 6 8 10]; % 6 8 10 = prw1 cohort
for i=1:length(tumorind)
    bpprob = [];
    bpprob = cell2mat(bpsall(1,tumorind(i),endpoint(tumorind(i))));
    bpprob(isnan(bpprob))=[];
    [N,~] = histcounts(bpprob,'Normalization','probability');
    bpprobs{i} = N;
end
for i=1:length(bpprobs)
    maxl(i) = length(bpprobs{i});
end
maxl = max(maxl);  
for i=1:length(bpprobs)
    bpprobs{i}(length(bpprobs{i})+1:maxl) = 0; %pad non representative numbers with zeros
end
bpprobs=cell2mat(bpprobs');
cdfs=zeros(size(bpprobs));
for j=1:length(bpprobs(:,1))
    cdfs(j,:) = cumsum(bpprobs(j,:));
end
meancdfs1=mean(cdfs,1);
%reset for next timepoint
clear bpprobs;
for i=1:length(tumorind)
    bpprob = [];
    bpprob = cell2mat(bpsall(1,tumorind(i),endpoint(tumorind(i))-2));
    bpprob(isnan(bpprob))=[];
    [N,~] = histcounts(bpprob,'Normalization','probability');
    bpprobs{i} = N;
end
for i=1:length(bpprobs)
    maxl(i) = length(bpprobs{i});
end
maxl = max(maxl);  
for i=1:length(bpprobs)
    bpprobs{i}(length(bpprobs{i})+1:maxl) = 0;
end
bpprobs=cell2mat(bpprobs');
cdfs=zeros(size(bpprobs));
for j=1:length(bpprobs(:,1))
    cdfs(j,:) = cumsum(bpprobs(j,:));
end
meancdfs2=mean(cdfs,1);
%reset for next timepoint
clear bpprobs;
for i=1:length(tumorind)
    bpprob = [];
    bpprob = cell2mat(bpsall(1,tumorind(i),endpoint(tumorind(i))-7));
    bpprob(isnan(bpprob))=[];
    [N,~] = histcounts(bpprob,'Normalization','probability');
    bpprobs{i} = N;
end
for i=1:length(bpprobs)
    maxl(i) = length(bpprobs{i});
end
maxl = max(maxl);  
for i=1:length(bpprobs)
    bpprobs{i}(length(bpprobs{i})+1:maxl) = 0;
end
bpprobs=cell2mat(bpprobs');
cdfs=zeros(size(bpprobs));
for j=1:length(bpprobs(:,1))
    cdfs(j,:) = cumsum(bpprobs(j,:));
end
meancdfs21=mean(cdfs,1);
%reset for next timepoint
clear bpprobs;
for i=1:length(tumorind)
    bpprob = [];
    bpprob = cell2mat(bpsall(1,tumorind(i),endpoint(tumorind(i))-15));
    bpprob(isnan(bpprob))=[];
    [N,~] = histcounts(bpprob,'Normalization','probability');
    bpprobs{i} = N;
end
for i=1:length(bpprobs)
    maxl(i) = length(bpprobs{i});
end
maxl = max(maxl);  
for i=1:length(bpprobs)
    bpprobs{i}(length(bpprobs{i})+1:maxl) = 0;
end
bpprobs=cell2mat(bpprobs');
cdfs=zeros(size(bpprobs));
for j=1:length(bpprobs(:,1))
    cdfs(j,:) = cumsum(bpprobs(j,:));
end
meancdfs3=mean(cdfs,1);
clear bpprobs;
for i=1:length(tumorind)
    bpprob = [];
    bpprob = cell2mat(bpsall(1,tumorind(i),endpoint(tumorind(i))-30));
    bpprob(isnan(bpprob))=[];
    [N,~] = histcounts(bpprob,'Normalization','probability');
    bpprobs{i} = N;
end
for i=1:length(bpprobs)
    maxl(i) = length(bpprobs{i});
end
maxl = max(maxl);  
for i=1:length(bpprobs)
    bpprobs{i}(length(bpprobs{i})+1:maxl) = 0;
end
bpprobs=cell2mat(bpprobs');
cdfs=zeros(size(bpprobs));
for j=1:length(bpprobs(:,1))
    cdfs(j,:) = cumsum(bpprobs(j,:));
end
meancdfs4=mean(cdfs,1);
clear bpprobs;
ctrlind = [4 5];
for i=1:length(ctrlind)
    bpprob = [];
    bpprob = cell2mat(bpsall(1,ctrlind(i),endpoint(ctrlind(i))));
    bpprob(isnan(bpprob))=[];
    [N,~] = histcounts(bpprob,'Normalization','probability');
    bpprobs{i} = N;
end
for i=1:length(bpprobs)
    maxl(i) = length(bpprobs{i});
end
maxl = max(maxl);  
for i=1:length(bpprobs)
    bpprobs{i}(length(bpprobs{i})+1:maxl) = 0;
end
bpprobs=cell2mat(bpprobs');
cdfs=zeros(size(bpprobs));
for j=1:length(bpprobs(:,1))
    cdfs(j,:) = cumsum(bpprobs(j,:));
end
meancdfs5=mean(cdfs,1); %control timepoint 1
clear bpprobs;
for i=1:length(ctrlind)
    bpprob = [];
    bpprob = cell2mat(bpsall(1,ctrlind(i),endpoint(ctrlind(i))-2));
    bpprob(isnan(bpprob))=[];
    [N,~] = histcounts(bpprob,'Normalization','probability');
    bpprobs{i} = N;
end
for i=1:length(bpprobs)
    maxl(i) = length(bpprobs{i});
end
maxl = max(maxl);  
for i=1:length(bpprobs)
    bpprobs{i}(length(bpprobs{i})+1:maxl) = 0;
end
bpprobs=cell2mat(bpprobs');
cdfs=zeros(size(bpprobs));
for j=1:length(bpprobs(:,1))
    cdfs(j,:) = cumsum(bpprobs(j,:));
end
meancdfs6=mean(cdfs,1); %control timepoint 2
clear bpprobs
for i=1:length(ctrlind)
    bpprob = [];
    bpprob = cell2mat(bpsall(1,ctrlind(i),endpoint(ctrlind(i))-15));
    bpprob(isnan(bpprob))=[];
    [N,~] = histcounts(bpprob,'Normalization','probability');
    bpprobs{i} = N;
end
for i=1:length(bpprobs)
    maxl(i) = length(bpprobs{i});
end
maxl = max(maxl);  
for i=1:length(bpprobs)
    bpprobs{i}(length(bpprobs{i})+1:maxl) = 0;
end
bpprobs=cell2mat(bpprobs');
cdfs=zeros(size(bpprobs));
for j=1:length(bpprobs(:,1))
    cdfs(j,:) = cumsum(bpprobs(j,:));
end
meancdfs7=mean(cdfs,1); %control timepoint 3
clear bpprobs
for i=1:length(ctrlind)
    bpprob = [];
    bpprob = cell2mat(bpsall(1,ctrlind(i),endpoint(ctrlind(i))-30));
    bpprob(isnan(bpprob))=[];
    [N,~] = histcounts(bpprob,'Normalization','probability');
    bpprobs{i} = N;
end
for i=1:length(bpprobs)
    maxl(i) = length(bpprobs{i});
end
maxl = max(maxl);  
for i=1:length(bpprobs)
    bpprobs{i}(length(bpprobs{i})+1:maxl) = 0;
end
bpprobs=cell2mat(bpprobs');
cdfs=zeros(size(bpprobs));
for j=1:length(bpprobs(:,1))
    cdfs(j,:) = cumsum(bpprobs(j,:));
end
meancdfs8=mean(cdfs,1); %control timepoint 4

fig2 = figure();
plot(meancdfs1,"Color",[0 0.5 0],"LineWidth",3); hold on
plot(meancdfs2,"Color",[0 0.7 0],"LineWidth",3); hold on
plot(meancdfs21,"Color",[0 0.8 0],"LineWidth",3); hold on
plot(meancdfs3,"Color",[0 0.9 0],"LineWidth",3); hold on
plot(meancdfs4,"Color",[0 1 0],"LineWidth",3); hold on
plot(meancdfs5,"Color",[0 0 0],"LineWidth",3); hold on
plot(meancdfs6,"Color",[0.5 0.5 0.5],"LineWidth",3); hold on
plot(meancdfs7,"Color",[0.7 0.7 0.7],"LineWidth",3); hold on
plot(meancdfs8,"Color",[0.9 0.9 0.9],"LineWidth",3); hold on
legend("last day", "last day-2", "last day-10","last day-15", "last day-30",...
    "last day", "last day-2", "last day-15", "last day-30");
xlabel("breaking point")
ylabel("probability")
title("cdf of breaking point over training and cancer progression")
text(25,0.2,'tumor group, n=6, control group n=2');
currfile = strcat(dst, '/', "cdf_over_time_prw1_prw2.fig");
saveas(fig2, currfile)
%%
fig3 = figure();
subplot(241);
% alignt one with the lowest endpoint
%mean across groups
tumor = [reshape(mbpall(1, 1, endpoint(1)-min(endpoint)+1:endpoint(1)), 1, []);...
    reshape(mbpall(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, []);...
    reshape(mbpall(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, [])];
%     reshape(mbpall(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []);...
%   reshape(mbpall(1, 8, endpoint(8)-min(endpoint)+1:endpoint(8)), 1, []);reshape(mbpall(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, [])];
tumormean = mean(tumor,1);
tumorstd = std(tumor, 1)/sqrt(length(tumor(:,1)));
ctrl = [reshape(mbpall(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []);...
    reshape(mbpall(1, 5, endpoint(5)-min(endpoint)+1:endpoint(5)), 1, [])];
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
    reshape(relaall(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, []);...
    reshape(relaall(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, [])];
%     reshape(relaall(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []);...
%   reshape(relaall(1, 8, endpoint(8)-min(endpoint)+1:endpoint(8)), 1, []);reshape(relaall(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, [])];
tumormean = mean(tumor,1);
tumorstd = std(tumor, 1)/sqrt(length(tumor(:,1)));
ctrl = [reshape(relaall(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []);...
    reshape(relaall(1, 5, endpoint(5)-min(endpoint)+1:endpoint(5)), 1, [])];
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
    reshape(rewaall(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, []); ...
    reshape(rewaall(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, []); ...
    reshape(rewaall(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, [])];
%     reshape(rewaall(1, 8, endpoint(8)-min(endpoint)+1:endpoint(8)), 1, []); ...
%     reshape(rewaall(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, [])];
tumormean = mean(tumor,1);
tumorstd = std(tumor, 1)/sqrt(length(tumor(:,1)));
ctrl = [reshape(rewaall(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []);...
    reshape(rewaall(1, 5, endpoint(5)-min(endpoint)+1:endpoint(5)), 1, [])];
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
    reshape(pokeall(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, []); ...
reshape(pokeall(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, [])];
%     reshape(pokeall(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []);...
%     reshape(pokeall(1, 8, endpoint(8)-min(endpoint)+1:endpoint(8)), 1, []);reshape(pokeall(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, [])];
tumormean = mean(tumor,1);
tumorstd = std(tumor, 1)/sqrt(length(tumor(:,1)));
ctrl = [reshape(pokeall(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []);...
    reshape(pokeall(1, 5, endpoint(5)-min(endpoint)+1:endpoint(5)), 1, [])];
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
    reshape(rxinits(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, []);...
    reshape(rxinits(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, [])];
%     reshape(rxinits(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []);...
%     reshape(rxinits(1, 8, endpoint(8)-min(endpoint)+1:endpoint(8)), 1, []); ...
%     reshape(rxinits(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, [])];
tumormean = mean(tumor,1);
tumorstd = std(tumor, 1)/sqrt(length(tumor(:,1)));
ctrl = [reshape(rxinits(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []);...
    reshape(rxinits(1, 5, endpoint(5)-min(endpoint)+1:endpoint(5)), 1, [])];
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
    reshape(rxsides(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, []);...
    reshape(rxsides(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, [])];
%     reshape(rxsides(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []);...
%     reshape(rxsides(1, 8, 1:endpoint(8)), 1, []);reshape(rxsides(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, [])];
tumormean = mean(tumor,1);
tumorstd = std(tumor, 1)/sqrt(length(tumor(:,1)));
ctrl = [reshape(rxsides(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []);...
    reshape(rxsides(1, 5, endpoint(5)-min(endpoint)+1:endpoint(5)), 1, [])];
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
    reshape(inactivesessions(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, []);...
    reshape(inactivesessions(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, [])];
%     reshape(inactivesessions(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []);...
%     reshape(inactivesessions(1, 8, endpoint(8)-min(endpoint)+1:endpoint(8)), 1, []);reshape(inactivesessions(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, [])];
tumormean = mean(tumor,1);
tumorstd = std(tumor, 1)/sqrt(length(tumor(:,1)));
ctrl = [reshape(inactivesessions(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []);...
    reshape(inactivesessions(1, 5, endpoint(5)-min(endpoint)+1:endpoint(5)), 1, [])];
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
text(0,0.5,'tumor group, n=3, control group n=2');
text(0,0.3, 'error bars represent standard error');
set(ax,'visible','off')
sgtitle('behavior variables for progressive ratio task')
currfile = strcat(dst, '/', "behavior_variables_summary_prw2.fig");
saveas(fig, currfile)
%%

fig = figure();
subplot(241);
% alignt one with the lowest endpoint
plot(reshape(mbpall(1, 1, endpoint(1)-min(endpoint)+1:endpoint(1)), 1, []), 'g', 'LineWidth',2), hold on
plot(reshape(mbpall(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, []), 'g', 'LineWidth',2), hold on
plot(reshape(mbpall(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(mbpall(1, 8, :), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(mbpall(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(mbpall(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(mbpall(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(mbpall(1, 7, endpoint(7)-min(endpoint)+1:endpoint(7)), 1, []), '--k', 'LineWidth',2), hold on
plot(reshape(mbpall(1, 9, endpoint(9)-min(endpoint)+1:endpoint(9)), 1, []), '--k', 'LineWidth',2), hold on
xticks([1:8:min(endpoint),32])
xticklabels(["back",-24,-16,-8,0])
ylabel('mean breaking point')
title('mean breaking point')
xlabel('days before endpoint')
legend('LHW2', 'NHW2', 'LHW1', 'NHW1', 'RHRHW1', ...
    'RHRHW2', 'RHW2', 'LHRHW1', 'RHW1', 'Location','northwest','NumColumns',2)
hold off

subplot(242);
% alignt one with the lowest endpoint
plot(reshape(relaall(1, 1, endpoint(1)-min(endpoint)+1:endpoint(1)), 1, []), 'g', 'LineWidth',2), hold on
plot(reshape(relaall(1, 3, endpoint(3)-min(endpoint)+1:endpoint(3)), 1, []), 'g', 'LineWidth',2), hold on
plot(reshape(relaall(1, 6, endpoint(6)-min(endpoint)+1:endpoint(6)), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(relaall(1, 8, :), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(relaall(1, 10, endpoint(10)-min(endpoint)+1:endpoint(10)), 1, []), '--g', 'LineWidth',2), hold on
plot(reshape(relaall(1, 4, endpoint(4)-min(endpoint)+1:endpoint(4)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(relaall(1, 2, endpoint(2)-min(endpoint)+1:endpoint(2)), 1, []), 'k', 'LineWidth',2), hold on
plot(reshape(relaall(1, 7, endpoint(7)-min(endpoint)+1:endpoint(7)), 1, []), '--k', 'LineWidth',2), hold on
plot(reshape(relaall(1, 9, endpoint(9)-min(endpoint)+1:endpoint(9)), 1, []), '--k', 'LineWidth',2), hold on
xticks([1:8:min(endpoint),32])
xticklabels(["back",-24,-16,-8,0])
ylabel('relation high/low')
xlabel('days before endpoint')
title('relation high/low rewards')
legend('LHW2', 'NHW2', 'LHW1', 'NHW1', 'RHRHW1', ...
    'RHRHW2', 'RHW2', 'LHRHW1', 'RHW1', 'Location','northwest','NumColumns',2)
hold off

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
title('trial activity for RHRHW2')
xlabel('trial number')
ylabel('CurRewardNum')
legend('last day', 'last day-1', 'last day-2', 'last day-20')
subplot(337)
plot(cell2mat(rewnums(1,5,endpoint(5))),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,5,endpoint(5)-1)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,5,endpoint(5)-2)),'LineWidth',2),hold on
plot(cell2mat(rewnums(1,5,endpoint(5)-20)),'k','LineWidth',2),hold on
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
