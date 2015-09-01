function [] = bs_limo_display_1stReduced(predictor,flagsOne,interaction)
subjStruct = flagsOne.subj;
dat = [];

if interaction
    FileName = sprintf('Interaction_effect_%i.mat',predictor)
else
    FileName = sprintf('Condition_effect_%i.mat',predictor)
end

for s = 1:length(subjStruct)
    % get clusterTimes
   cd(subjStruct(s).folder);
   load(subjStruct(s).LIMO);
   alpha = 0.05;
   MCC = 2;
   h0File = sprintf(['%s/H0/mcc_%i_%s'],subjStruct(s).folder,MCC,FileName);
   fprintf('assuming clusterTest MCC=2 \n')
   if exist(h0File,'file')
        load(h0File)
    
   else
    
        [M, mask, mytitle] = limo_stat_values(1,FileName,alpha,MCC,LIMO,'use empirical p values',[]);
        if isnan(M)
            M = nan(size(mask));
        end
        
        
        if MCC>1 && ~isempty(mask)
            save(h0File,'mask','M','mytitle')
        end
        
    end
    load(sprintf('%s/%s',subjStruct(s).folder,FileName))
    if interaction
        pMatrix =  squeeze(Interaction_effect(:,:,2));
    else
    pMatrix =  squeeze(Condition_effect(:,:,2));
    end
    
    sigTimes = any(isnan(M));
    
    
    sigVals = nan(1,400);
    for t = 1:400
        if sigTimes(t)
        tmpVal = min(M(:,t));
        if ~isempty(tmpVal); 
            sigVals(t) = tmpVal;    
        end
        end
    end
    dat.allSigTimes(s,:) = sigTimes;
    dat.allSigVals(s,:) = sigVals;
    
end


%%
timevect = linspace(LIMO.data.start*1000,LIMO.data.end*1000,size(dat.allSigTimes,2));
dat.allSigVals(isnan(dat.allSigVals)) = 0.05+eps;
figure
imagesc(timevect,1:size(dat.allSigTimes,1),dat.allSigVals)
% set('YTick',1:size(dat.allSigTimes,1))
colormap([cbrewer('seq','RdPu',100,[],1);0.7 0.7 0.7]);
caxis([0 0.05])
colorbar
% title(titleString)
title(sprintf('%s    -    %s',flagsOne.name,LIMO.design.XDesc{predictor}))