function [] = bs_limo_topoplots(flagIdx,param)

flags = be_check_folderstruct('bsLimo');
p = flags(flagIdx);

if length(param) == 1
   param(2) = param(1)*2;
   param(1) = param(2)-1;
end

fileName = sprintf('paired_samples_ttest_parameter_%i%i.mat',param(1),param(2));
cd(p.group.folder);
load(p.group.LIMO);
load(fileName);

%% Load Cluster Mass
MCC = 5;
if exist(sprintf('H0/mcc_%i_%s',MCC,fileName),'file')
    load(sprintf('H0/mcc_%i_%s',MCC,fileName))
end
if ~exist('mask','var') || isempty(mask)
    
    [M, mask, mytitle] = limo_stat_values(1,fileName,p,MCC,LIMO, 'use theoretical p values',[]);
    if isnan(M)
        M = nan(size(mask));
    end
    
    
    if MCC>1 && ~isempty(mask)
        save(sprintf('H0/mcc_%i_%s',MCC,fileName),'mask','M','mytitle')
    end
    
end

%% Calculate data to topoplot
% topoData = squeeze(paired_samples(:,:,4));
% pMatrix =  squeeze(paired_samples(:,:,5));
rawDataS = bs_limoGatherData('paths',LIMO.data.fullPath,'predictor',param);
rawData = mean(rawDataS,4);
timevec = round((LIMO.data.start:1/LIMO.data.sampling_rate:LIMO.data.end)*1000);
topoDiff = rawData(:,:,1) - rawData(:,:,2);
figure
scale = max([abs(min(topoDiff(:))) abs(max(topoDiff(:)))]);
topoTimes = [-300:25:400];

jisubplot(5,8,0,'landscape',[0.1 0.2])
for k = 1:length(topoTimes)-1
    nextplot('byrow')
    fromT = find(timevec<=topoTimes(k),1,'last');
    toT = find(timevec<topoTimes(k+1),1,'last');
    tmpDat= mean(topoDiff(:,fromT:toT),2);
    topoplot(tmpDat,LIMO.data.chanlocs,'maplimits',[-scale scale]);
    title(sprintf('%.0f <= %.0f',timevec(fromT),timevec(toT)))
    
    
end

cbar