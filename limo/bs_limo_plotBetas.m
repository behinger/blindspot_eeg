LIMO.dir = '/net/store/nbp/EEG/blind_spot/data/limo01-07';
[~,~,fullPath] = bs_limo_getFilepaths('filename','Betas.mat','subjects',1:10,'path',LIMO.dir);
data =bs_limoGatherData('paths',fullPath);

eegTime = linspace(-300,498,400);
%%
figure

elec = 36;
subj = 1;
leg = [];
colorSel = [1 0 0; 0 1 1; 0 0 1; 0 1 0];%jet(4);

for k = 1:2:8
    leg(ceil(k/2)) = plot(eegTime,squeeze(data(elec,:,k,subj)),'Color',colorSel(ceil(k/2),:)); hold on
    plot(eegTime,squeeze(data(elec,:,k+1,subj)),'--','Color',colorSel(ceil(k/2),:)); hold all
end
leg(end+1) = plot(eegTime,squeeze(data(elec,:,end,subj)),'Color',[0 0 0]); hold on
legend(leg,'eye','bs','temp','spat','const')
% plot(squeeze(Betas(36,:,k+1)),'-.'),hold all

hline(0,'k')
axis([ -50 400 -3 3])

%% mean Betas over subjects


elec = 36;
legendList = {'eye','bs','temp','spat','const'};
colorSel = winter(5);
colorSubj = jet(10);
axisLim = [ -150 400 -4 6];
for k = [1:2:8 9]
    figure
    subplot(3,1,[1]), % 4 7
    %     limo_trimci
    dataCI = bootci(100,@mean,squeeze(data(elec,:,k,:))')';
    %     dataCI = squeeze(dataCI(elec,:,:));
    CIs = abs(bsxfun(@minus,dataCI(:,:),mean(squeeze(data(elec,:,k,:)),2)));
    leg(ceil(k/2)) = boundedline(eegTime,dataCI(:,2),CIs(:,[1 2]),'cmap',colorSel(ceil(k/2),:),'alpha');
    if k ~=9
        %     dataCI = limo_trimmed_mean(squeeze(data(:,:,k+1,:)),1,0.05);
        %     dataCI = squeeze(dataCI(elec,:,:));
        %     CIs = abs(bsxfun(@minus,dataCI(:,[1 3]),dataCI(:,2)));
        dataCI = bootci(100,@mean,squeeze(data(elec,:,k+1,:))')';
        CIs = abs(bsxfun(@minus,dataCI(:,:),mean(squeeze(data(elec,:,k+1,:)),2)));
        boundedline(eegTime,dataCI(:,2),CIs(:,[1 2]),'--','cmap',colorSel(ceil(k/2),:)*0.5,'alpha');
    end
    axis(axisLim)
    vline(0),hline(0,'k--')
    title(sprintf('Elec %i, pred: %s',elec,legendList{ceil(k/2)}))
    set(gca,'XTick',[])

    subplot(3,1,2), % 7
    for subj = 1:10
       h1 =  patchline(eegTime,squeeze(data(elec,:,k,subj)),'EdgeAlpha',0.5);
        set(h1,'EdgeColor',colorSubj(subj,:));

        if k ~=9
            %         h2(subj) = patchline(eegTime,squeeze(data(elec,:,k+1,subj)),'EdgeAlpha',0.1);
            %         set(h2,'EdgeColor',colorSel(ceil(k/2),:)*0.5);
        end
    end
    axis(axisLim)
    vline(0),hline(0,'k--')
    set(gca,'XTick',[])

    subplot(3,1,3), % 7
    for subj = 1:10
%         patchline(eegTime,squeeze(data(elec,:,k,subj)),'EdgeAlpha',0.1)
        if k ~=9
                    h2 = patchline(eegTime,squeeze(data(elec,:,k+1,subj)),'EdgeAlpha',0.5);
                    set(h2,'EdgeColor',colorSubj(subj,:));
        end
    end
    axis(axisLim)
    vline(0),hline(0,'k--')
    
    be_print('file',sprintf('/work/behinger/Dropbox/Masterarbeit/Plots/eeg/1stLevel/%s',legendList{ceil(k/2)}),'prepare','yes')
end

% plot(squeeze(Betas(36,:,k+1)),'-.'),hold all

% hline(0,'k')
%%
figure

selPred = 3;
legendList = {'eye','bs','temp','spat','const'};



timePlt = [-100:20:400];
colorSel = jet(4);
for k = 1:length(timePlt)
    subplot(5,6,k)
    dataPlot = squeeze(nanmean(nanmean(data(:,eegTime==timePlt(k),selPred,:),4),2));
    dataPlot2 = squeeze(nanmean(nanmean(data(:,eegTime==timePlt(k),selPred+1,:),4),2));
    topoplot(dataPlot-dataPlot2,LIMO.data.chanlocs,'maplimits',[-2 2],'electrodes','off');
    title(sprintf('timePlt:%i',timePlt(k)))
    
    
end
title(sprintf('Diff Pred %s',legendList{ceil(selPred/2)}))


% plot(squeeze(Betas(36,:,k+1)),'-.'),hold all
%%
figure

selPred = 3;
legendList = {'eye','eye','bsOff','bsOn','temp','temp','spatOff','spatOn','const'};
timePlt = [-100:20:400];
colorSel = jet(4);
for k = 1:length(timePlt)
    subplot(5,6,k)
    dataPlot = squeeze(nanmean(nanmean(data(:,eegTime==timePlt(k),selPred,:),4),2));
    topoplot(dataPlot,LIMO.data.chanlocs,'maplimits',[-2 2],'electrodes','off');
    title(sprintf('timePlt:%i',timePlt(k)))
    
    
end
title(sprintf('Pred %s',legendList{selPred}))% plot(squeeze(Betas(36,:,k+1)),'-.'),hold all
%% Mean R2 topoplot
[~,~,fullPath] = bs_limo_getFilepaths('filename','R2.mat','subjects',1:10,'path',LIMO.dir);
data =bs_limoGatherData('paths',fullPath);

timePlt = [-100:20:400];
subjList = {1 2 3 4 5 6 7 8 9 10 1:10}
for l = length(subjList)
    figure
    subj = subjList{l}
    for k = 1:length(timePlt)
        subplot(5,6,k)
        dataPlot = squeeze(nanmean(data(:,eegTime==timePlt(k),subj),3));
        topoplot(dataPlot,LIMO.data.chanlocs,'maplimits',[0 0.05],'electrodes','off');
        title(sprintf('timePlt:%i',timePlt(k)))
        
    end
    cbar
    title(sprintf('R2 Subj %i',subj))% plot(squeeze(Betas(36,:,k+1)),'-.'),hold all
end