function [] = bs_limo_plotEEGinteraction(flags,flagIdx,pred,method)


cfg = [];
if nargin < 4
    method = '';
end
cfg.norm = strcmp(method,'normalized');
cfg.flagIdx = flagIdx;
cfg.pred = pred;
FileName = sprintf('one_sample_ttest_parameter_%i.mat',cfg.pred);
cd(flags(cfg.flagIdx).group.folder)
load ([FileName]);

load(flags(cfg.flagIdx).group.LIMO)
betas = squeeze(one_sample(:,:,4));
% pMatrix =  squeeze(one_sample(:,:,5));
choice = 'use theoretical p values';
[~,LI] = bs_limo_plotDesignmatrix(LIMO);
cfg.xDesc = bs_limo_designMat_betterDesc(LI.design.XDesc);
cfg.timevec = round((LIMO.data.start:1/LIMO.data.sampling_rate:LIMO.data.end)*1000);
[M, mask, mytitle] = limo_stat_values(1,FileName,0.05,5,LIMO,choice,[]);
%% Four Predictors of interest:
startI = sum(LI.design.nb_conditions)+sum(LI.design.nb_interactions(1:cfg.pred-100-1));
endI =  LI.design.nb_interactions(cfg.pred-100);
cfg.interactPreds = startI+1:startI+endI;
cfg.mainPred = [];
for k = cfg.interactPreds
    diffSet = [];
for mE1 = 1:sum(LI.design.nb_conditions) %mainEffect
    % Two Way
    for mE2 = mE1+1:sum(LI.design.nb_conditions)
        if all(LI.design.X(:,k) == LI.design.X(:,mE1).*LI.design.X(:,mE2))
            diffSet = [diffSet {[mE1 mE2]}];
        end
    end
end
cfg.mainPred = [cfg.mainPred {diffSet}];
end

fprintf('loading started \n')
clear rawDataOrg;
load([flags(cfg.flagIdx).group.folder '/YrAll.mat'])
rawDataOrg = YrAll(:,:,cfg.interactPreds);
% [~,pairs] = bs_limo_designMat_removeDouble(LI);
% rawDataOrg = bs_limoGatherData('paths',LIMO.data.fullPath,'predictor',cfg.interactPreds);
% rawDataOrg = bs_limoGatherData('paths',cellfun(@(x)[ x '/Yr.mat'],LIMO.data.data_dir,'UniformOutput',0),'predictor',cfg.interactPreds);


% doublePreds = cfg.interactPreds(ismember(cfg.interactPreds,pairs));
% for w= 1:length(doublePreds)
%     for x = 1:length(pairs)
%         if any(doublePreds(w) == pairs(x,:))
%             newPred = pairs(x,mod(find(doublePreds(w) == pairs(x,:)),2)+1); % find and choose the other one
%             newRawDat = bs_limoGatherData('paths',LIMO.data.fullPath,'predictor',newPred);
%             newRawDat = permute(newRawDat,[1 2 4 3]); % just a trick to get the 3 dimension to be 1 for adding in the next line
%             rawDataOrg(:,:,doublePreds(w)==cfg.interactPreds,:) = rawDataOrg(:,:,doublePreds(w)==cfg.interactPreds,:)+newRawDat;
%         end
%     end
% end

fprintf('loading finished \n')
%%



figure

topoTimes = linspace(-100,450,16);

jisubplot(4 ,15,0,'landscape',[0.1 0.2])
predList = {[1 2] [3 4] [1 3] [2 4]};
% rawDatOrg(:,:,:,:) = rawDatO
for l = 1:4
    %find out which predictor is constant
    tmpMainPred = [cfg.mainPred{predList{l}}];
    tmpMainPred = [tmpMainPred{:}];
    [C,ia,ic] = unique(tmpMainPred);
    constPred = tmpMainPred(setdiff([1:4],ia));
    if cfg.norm
    topoDiff = mean(rawDataOrg(:,:,l,:),4) -mean(rawDataOrg(:,:,1,:),4)  ;
    else
    topoDiff = mean(rawDataOrg(:,:,predList{l}(2),:),4) -mean(rawDataOrg(:,:,predList{l}(1),:),4)  ;
    end
%     topoDiff = topoDiff - (YrAll
    scale = max([abs(min(topoDiff(:))) abs(max(topoDiff(:)))]);
    for k = 1:length(topoTimes)-1
        if k == 1; 
            nextplot('newrow');
            if cfg.norm
                title(sprintf('%s',[cfg.xDesc{[cfg.mainPred{l}{:}]}]))
            else
                title(sprintf('%s',cfg.xDesc{constPred}));
            end
            elec = 'ptslabels';
            
        else
            elec = 'on';
            nextplot('byrow');
            title(sprintf('%.0f <= %.0f',cfg.timevec(fromT),cfg.timevec(toT)))
        end
        fromT = find(cfg.timevec<=topoTimes(k),1,'last');
        toT = find(cfg.timevec<topoTimes(k+1),1,'last');
        tmpDat= mean(topoDiff(:,fromT:toT),2);
        topoplot(tmpDat,LIMO.data.chanlocs,'maplimits',[-scale scale],'nosedir','+Y','electrodes',elec);
        
    end
    
end

% be_tool_axesSameCaxis([-0.5 0.5])
be_tool_axesSameCaxis










% %%
% rawDat = rawDataOrg;
% clear interActValPos
% clear interActValNeg
% clear interActVal
% betaVals = repmat(betas,[1,1,15]);
% maskVP = repmat(mask,[1,1,15]);
% [~,c] = find(M(:,:)<=(0.05+eps));
% sigArea = unique(c)';
% if length(sigArea)>1
%     diffList = [diff([sigArea(1)-1 sigArea])~=1];
%     for k = 1:sum(diffList)+1
%         tmpIdx = sigArea(1);
%         if diffList(1) == 1
%             tmpEnd = tmpIdx;
%         else
%             tmpEnd = sigArea(find(diffList,1,'first')-1);
%         end
%         if isempty(tmpEnd)
%             tmpEnd = sigArea(end);
%         end
%         sigPart = [tmpIdx tmpEnd];
%         sigPart = timevect(sigPart)
%         for l = 1:4
%             rawDataMask = bsxfun(@times,mask,squeeze(rawData(l,:,:,:)));
%             rawDataMask(maskVP==0) = nan;
%             interActVal(k,l,:,:) = squeeze(nanmean((rawDataMask(:,tmpIdx:tmpEnd,:)),2));
%             
%             tmpRaw = rawDataMask(:,tmpIdx:tmpEnd,:);
%             tmpRaw(betaVals(:,tmpIdx:tmpEnd,:)<0) = nan;
%             interActValPos(k,l,:,:) = squeeze(nanmean((tmpRaw),2));
%             tmpRaw = rawDataMask(:,tmpIdx:tmpEnd,:);
%             tmpRaw(betaVals(:,tmpIdx:tmpEnd,:)>0) = nan;
%             interActValNeg(k,l,:,:) = squeeze(nanmean((tmpRaw),2));
%         end
%         
%         
%         if sum(diffList) ~= 0
%             diffList(find(sigArea==tmpIdx):find(sigArea==tmpEnd)+1) = [];
%             sigArea(find(sigArea==tmpIdx):find(sigArea==tmpEnd)+1) = [];
%         end
%     end
% end