function bs_limo_display_resultsV4(pred,PathName,sigElec,varargin)

g = be_inputcheck(varargin,...
    {'recalcMCC','boolean',[],0;
    'contour','string',{'yes','no'},'yes';
    'MCC','integer',[1:5],2;
    'alpha','real',[0 1],0.05;
    'interaction','string',{'yes','no'},'no';
    'style','string','','';
    'caxis','integer',[],[];
    'plotAxis','integer',[],[];
    'selChan','integer',[],[];
    'topoTime','integer',[],[]});
if ischar(g)
    error(g)
end

cd(PathName)
load('LIMO.mat')
[~,LI] = bs_limo_plotDesignmatrix(LIMO);
[~,a] = fileparts(LIMO.dir);
predOrg = pred;
if strcmp(g.interaction,'yes')
    FileName = sprintf('one_sample_ttest_parameter_%i.mat',pred);
    g.title = [a '--- Interaction:' num2str(pred) ];
else
    g.title = [a '---' LI.design.XDesc{pred} ];
    if length(pred) == 1
        pred(2) = pred(1)*2;
        pred(1) = pred(2)-1;
    end
    FileName = sprintf('paired_samples_ttest_parameter_%i%i.mat',pred);
end
load (FileName);

%% Define Predictors
if ischar(g)
    error(g)
end
choice = 'use theoretical p values';
LI.design.XDescAll = bs_limo_designMat_betterDesc(LI.design.XDesc);
[predName] = bs_limo_designmat_shortPredToLongPred(LI,predOrg);


if strcmp(g.interaction,'yes')
    l = pred-100;
    if l < 0
        l = 1;
    end
    startI = sum(LI.design.nb_conditions)+sum(LI.design.nb_interactions(1:l-1));
    endI =  LI.design.nb_interactions(l);
    preds = startI+1:startI+endI;
elseif length(FileName) == 39
    preds = [ str2num(FileName(32:end-6)) str2num(FileName(34:end-4))];
else
    
    preds = [str2num(FileName(end-5:end-5)) str2num(FileName(end-4:end-4))];
end


%% Load the Data

if strncmp(FileName,'one_sample',10)
    toplot = squeeze(one_sample(:,:,4));
%     pMatrix =  squeeze(one_sample(:,:,5));
    
elseif strncmp(FileName,'paired',6)
    toplot = squeeze(paired_samples(:,:,4));
%     pMatrix =  squeeze(paired_samples(:,:,5));
    assignin('base','T_values',toplot)
    
end

% We need to check whether some of the predictors are actually coded
% double, then we identify them and just sum them as everything in a GLM is
% linear anyway.
[~,pairs] = bs_limo_designMat_removeDouble(LI);
rawDataS = bs_limoGatherData('paths',LIMO.data.fullPath,'predictor',preds);

doublePreds = preds(ismember(preds,pairs));
for w= 1:length(doublePreds)
    for x = 1:length(pairs)
        if any(doublePreds(w) == pairs(x,:))
            newPred = pairs(x,mod(find(doublePreds(w) == pairs(x,:)),2)+1); % find and choose the other one
            newRawDat = bs_limoGatherData('paths',LIMO.data.fullPath,'predictor',newPred);
            newRawDat = permute(newRawDat,[1 2 4 3]); % just a trick to get the 3 dimension to be 1 for adding in the next line
            rawDataS(:,:,doublePreds(w)==preds,:) = rawDataS(:,:,doublePreds(w)==preds,:)+newRawDat;
        end
    end
end

% Add them up correctly
if  strcmp(g.interaction,'yes')
    if length(preds) == 4
        rawData = rawDataS(:,:,1,:) - rawDataS(:,:,2,:)- rawDataS(:,:,3,:) + rawDataS(:,:,4,:);
    elseif length(preds) == 8
        rawData = (rawDataS(:,:,1,:) - rawDataS(:,:,2,:)- rawDataS(:,:,3,:) + rawDataS(:,:,4,:))-(rawDataS(:,:,5,:) - rawDataS(:,:,6,:)- rawDataS(:,:,7,:) + rawDataS(:,:,8,:));
    elseif length(preds) == 16
        rawData = ((rawDataS(:,:,1,:) - rawDataS(:,:,2,:)- rawDataS(:,:,3,:) + rawDataS(:,:,4,:))-(rawDataS(:,:,5,:) - rawDataS(:,:,6,:)- rawDataS(:,:,7,:) + rawDataS(:,:,8,:))) - ...
            ((rawDataS(:,:,9,:) - rawDataS(:,:,10,:)- rawDataS(:,:,11,:) + rawDataS(:,:,12,:))-(rawDataS(:,:,13,:) - rawDataS(:,:,14,:)- rawDataS(:,:,15,:) + rawDataS(:,:,16,:)));
    end
else
    rawData = rawDataS(:,:,1,:) - rawDataS(:,:,2,:);
end
rawData = -rawData; % because we put 1 as off and 2 as on, I think it is more intuitive the other way round...

timevect = linspace(LIMO.data.start*1000,LIMO.data.end*1000,size(toplot,2));

%% Load the statistics
if g.MCC == 5
    [mask M] = bs_limo_tfceCalc(g.MCC,FileName,g.recalcMCC,LIMO,toplot);
elseif g.MCC == 2
    if exist(sprintf('H0/mcc_%i_%s',g.MCC,FileName),'file') && ~g.recalcMCC
        load(sprintf('H0/mcc_%i_%s',g.MCC,FileName))
    end
end


if ~exist('mask','var') || isempty(mask)
    Type = 1;
    [M, mask, mytitle] = limo_stat_values(Type,FileName,g.alpha,g.MCC,LIMO,choice,[]);
    if isnan(M)
        M = nan(size(mask));
    end
    if g.MCC>1 && ~isempty(mask) && g.MCC ~= 6
        save(sprintf('H0/mcc_%i_%s',g.MCC,FileName),'mask','M','mytitle')
    end
    
end

if g.MCC == 2
    M(isnan(M)) = 1;
end
%%

% do the figure
figure; set(gcf,'Color','w');
tmp = get(gcf,'Position');
set(gcf,'Position',[tmp(1:2)   1301         408])
bs_preparePlosBio(struct('column',2));
s.line = 4;
s.topo = 2;
jisubplot(s.line+1+2*s.topo,14,0,'',[0.05 0.05])



%% Define Colorbars
cMap.p = ([cbrewer('seq','RdPu',100,[],1)]);
cMap.b = cbrewer('div','RdYlBu',size(jet,1),[],1);%
cMap.lines = cbrewer('qual','Set1',8);%
%     freezeColors
cMap.cohenD = cbrewer('div','BrBG',100);%
    dataToPlot = winMean(rawData(:,:,:),3);

    [cluster, num] = bs_limo_ft_findcluster_posNeg(dataToPlot,M<=0.05,LIMO.data.neighbouring_matrix);

if ~strcmp(g.style,'topoPlotOnly')
    
    %% Plot the ERP at the siginficant electrods
    
    nextplot('newRow','size',[14 s.line])
    
    h = plot(timevect,dataToPlot','Color',[0.6 0.6 0.6]);
    
    if ~isempty(g.selChan)
        set(h(g.selChan),'Color',[1 0 0])
    end
    hold all;
    
    %% Mark all alpha < 0.05
    for e = 1:64
        [~,sigArea] = find(M(e,:)<=(0.05+eps));
        if length(sigArea)>1
            tmp = [diff([sigArea(1)-1 sigArea])~=1];
            for k = 1:sum(tmp)+1
                if isempty(tmp) %|| length(tmp) == 1
                    continue
                end
                
                tmpIdx = sigArea(1);
                if tmp(1) == 1
                    tmpEnd = tmpIdx;
                else
                    tmpEnd = sigArea(find(tmp,1,'first')-1);
                end
                if isempty(tmpEnd)
                    tmpEnd = sigArea(end);
                end
                sigPart = [tmpIdx tmpEnd];
                if sum(tmp) ~= 0
                    if length(tmp) == 1
                    tmp(sigArea==tmpIdx) = [];
                    sigArea(sigArea==tmpIdx) = [];
                    else
                    tmp(find(sigArea==tmpIdx):find(sigArea==tmpEnd)+1) = [];
                    sigArea(find(sigArea==tmpIdx):find(sigArea==tmpEnd)+1) = [];
                    end
                end

                sigPart(sigPart>size(rawData,2)) = size(rawData,2);
                if (sigPart(2) - sigPart(1)) < 5
                    
                    continue
                end
                plot(timevect(sigPart(1):sigPart(2)),dataToPlot(e,sigPart(1):sigPart(2)),'LineWidth',1,'Color',[0.2157 0.4941 0.7216])
            end
        end
    end
    
    %% Plot the Bonferroni corrected ones
    if g.alpha <0.05 %do it only if you need to
        for e = 1:64
            [~,sigArea] = find(M(e,:)<=(g.alpha+eps));
            if length(sigArea)>1
                tmp = [diff([sigArea(1)-1 sigArea])~=1];
                for k = 1:sum(tmp)+1
                    if isempty(tmp) || length(tmp) == 1
                        continue
                    end
                    
                    tmpIdx = sigArea(1);
                    if tmp(1) == 1
                        tmpEnd = tmpIdx;
                    else
                        tmpEnd = sigArea(find(tmp,1,'first')-1);
                    end
                    if isempty(tmpEnd)
                        tmpEnd = sigArea(end);
                    end
                    sigPart = [tmpIdx tmpEnd];
                    if sum(tmp) ~= 0
                        tmp(find(sigArea==tmpIdx):find(sigArea==tmpEnd)+1) = [];
                        sigArea(find(sigArea==tmpIdx):find(sigArea==tmpEnd)+1) = [];
                    end
                    
                    sigPart(sigPart>size(rawData,2)) = size(rawData,2);
                    % We don't want to skip them, because we are interested
                    % where exactly they are
%                       if (sigPart(2) - sigPart(1)) < 5 
%                     continue
%                 end
                    plot(timevect(sigPart(1):sigPart(2)),dataToPlot(e,sigPart(1):sigPart(2)),'LineWidth',1,'Color',[0 0 0])
                end
            end
        end
    end
    
    
    %% Write down where significant effects start
    
    fprintf('\n\n %s \n Beta: %s, Alpha:%.4f \n',PathName,[predName],g.alpha)
    
%     [cluster, num] = limo_ft_findcluster(mask,  LIMO.data.neighbouring_matrix);
tmpSort = [];
    for k = 1:num
        tmpSort(k) = find(cluster == k,1,'first');
    end
    [~, num] = sort(tmpSort);
    for k = num
        
        sigPart = [ timevect(find(any(cluster==k,1),1,'first'))  timevect(find(any(cluster==k,1),1,'last')) ];
        medianP = median(M(cluster==k));
        minP = min(M(cluster==k));
        
        if diff(sigPart) < 10
            cluster(cluster == k) = 0;
            continue
        end
        if sum(sum(dataToPlot(cluster==k))) > 0
            posNeg = 'positive';
        else
            posNeg = 'negative';
        end
        fprintf('%s effect from %.0fms to %.0fms (median-p =%.3f, min-p = %.3f) \n',posNeg,sigPart,medianP,minP);
%         fprintf('%s effect from %.0fms to %.0fms (all p < 0.05 TFCE corrected, median-p =%.3f) \n',posNeg,sigPart,medianP);
        
    end
    %% General Plot Stuff
    ax = gca;
    set(ax,'Box','off')
    set(ax,'XTick',[-300:100:400],'YAxisLocation','right');
    if isempty(g.topoTime)
        
        set(ax,'XLim',[-275 425])
    else
        set(ax,'XLim',[g.topoTime(1) g.topoTime(2)])
    end
    
    
    tmpLabel = num2cell([-300:100:400]);
    % for k = length(tmpLabel):-1:1
    %    tmpLabel = [tmpLabel(1:k-1) repmat({''},1,3) tmpLabel(k:end)] ;
    % end
    set(ax,'XTickLabel',tmpLabel)
    if isempty(g.plotAxis)
        yLimCustom(1) = min(dataToPlot(:));
        yLimCustom(2) = max(dataToPlot(:));
        set(ax,'YLim',yLimCustom*1.1)
    else
        set(ax,'YLim',g.plotAxis)
    end
    
    xlabel('Time [ms]')
    
    ylabel(ax,'Betas [mV]')
    
    % Add Minor Ticks
    xg = get(ax,'XLim');
    xg = linspace(xg(1),xg(2),29);
    yg = get(ax,'YLim');
    yg = yg(1)+[0 0.025*(yg(2)-yg(1))];
    xx = reshape([xg;xg;NaN(1,length(xg))],1,length(xg)*3);
    yy = repmat([yg NaN],1,length(xg));
    h_minorgrid = plot(ax,xx,yy,'k');
    
    
    
    if strcmp(g.interaction,'yes')
        startIdx = length(LI.design.nb_conditions);
    else
    end
    
    hline(0,'k')
    vline(0,'k')
    
    title(sprintf('Beta: %s, Alpha:%.4f',predName,g.alpha))
    drawnow
    
else
  
    tmpSort = [];
    for k = 1:num
        tmpSort(k) = find(cluster == k,1,'first');
    end
    [~, num] = sort(tmpSort);
    for k = num
        
        sigPart = [ timevect(find(any(cluster==k,1),1,'first'))  timevect(find(any(cluster==k,1),1,'last')) ];
        medianP = median(M(cluster==k));
        minP = min(M(cluster==k));
        
        if diff(sigPart) < 10
            cluster(cluster == k) = 0;
            continue
        end
        if sum(sum(dataToPlot(cluster==k))) > 0
            posNeg = 'positive';
        else
            posNeg = 'negative';
        end
        fprintf('%s effect from %.0fms to %.0fms (median-p =%.3f, min-p = %.3f) \n',posNeg,sigPart,medianP,minP);
        %         fprintf('%s effect from %.0fms to %.0fms (all p < 0.05 TFCE corrected, median-p =%.3f) \n',posNeg,sigPart,medianP);
        
    end
end


%% Plot the Topoplots

% nextplot('size',[10 20])
topoDiff = squeeze(mean(rawData(:,:,:,:),4));
topoEffect = squeeze(mean(rawData(:,:,:,:),4))./squeeze(std(rawData(:,:,:,:),[],4));
if isempty(g.caxis)
    scale = max([abs(min(topoDiff(:))) abs(max(topoDiff(:)))]);
    scale = prctile(topoDiff(:),[1 99]);
    scale = max(abs(scale));
else
    scale = abs(max(g.caxis));
    
end
if ~isempty(g.topoTime)
    topoTimes = linspace(g.topoTime(1),g.topoTime(2),15);
else
    topoTimes = linspace(-275,425,15);
end


nextplot('newrow','size',[1,1])
set(gca,'Visible','off')
nextplot('newrow','size',[1,s.topo])
for k = 1:length(topoTimes)-1
    fromT = find(timevect<=topoTimes(k),1,'last');
    toT = find(timevect<topoTimes(k+1),1,'last');
    tmpDat= mean(topoDiff(:,fromT:toT),2);
    
    sigElecs = any(cluster(:,fromT:toT)~=0,2);
    
%     sigElecs = any((M(:,fromT:toT)<=0.05)');
    hcolor = 'k';
    if k == 1 && ~isempty(g.selChan)
        sigElecs(g.selChan) =1;
    end
    emarker2 = {[find(sigElecs)] '.','k',10,1};
    style = 'both';
    if strcmp(g.style,'onlyCurves')
        style = 'contour';
        
    elseif strcmp(g.style,'onlyMaps')
        style = 'map';
        hcolor = 'none';
        emarker2 = {[]};
        
    end
    
    topoplot(tmpDat,LIMO.data.chanlocs,'style',style,'maplimits',[-scale scale],...
        'electrodes','off','emarker2',emarker2,'whitebk','on','numcontour',[linspace(-1.1*scale,1.1*scale,8)],...
        'nosedir','+Y','headrad',0.8,'hcolor',hcolor,'plotrad',0.9);
    
    
    colormap(cMap.b),
    if ~strcmp(g.style,'onlyCurves')
        freezeColors;
    end
    
    if ~strcmp(g.style,'topoPlotOnly')
        if k == length(topoTimes)-1
            cLimits = caxis;
            c = cbar('vert',0,cLimits,3);
            newPos = get(c,'Position');
            newPos(1) = newPos(1)-0.018;
            yTicksCbar = [-floorTo(cLimits(2),-2) 0 floorTo(cLimits(2),-2)];
            set(c,'Position',newPos,'YTick',yTicksCbar,'YTickLabel',yTicksCbar)
            
            title(c,'[mV]')
            freezeColors(c)
        end
        
    end
    
    %             %% Effectsize
    %     nextplot('bycol')
    %     effectSize = mean(topoEffect(:,fromT:toT),2);
    %     %     pvalDat(pvalDat>0.05) = nan%0.05+eps;
    %     topoplot(effectSize,LIMO.data.chanlocs,'style',style,'maplimits',[-4 4],...
    %         'electrodes','off','emarker2',emarker2,...
    %         'nosedir','+Y','headrad',0.8,'hcolor',hcolor,'plotrad',0.9);
    % %     caxis([-3 log10(0.5)+eps])
    %     if ~strcmp(g.style,'onlyCurves')
    %
    %         colormap([cMap.cohenD]),
    %         freezeColors;
    %     end
    %     if ~strcmp(g.style,'topoPlotOnly') && k == length(topoTimes)-1
    %
    %         c = cbar;
    %         cbfreeze(c);
    %     end
    %% Pvalue
    nextplot('bycol','size',[1 s.topo])
    %     pvalDat = mean(M(:,fromT:toT),2);
    pvalDat = min(M(:,fromT:toT),[],2);
    
    if any(pvalDat(:)<=0.05)
    pvalDat = log10(pvalDat);
    pvalDat(pvalDat==-Inf) = -4;
%     tmp = pvalDat;
%     tmp(pvalDat>log10(0.05)) = log10(0.05);
%     num_contur = min(length(unique(tmp)),6);
    %     pvalDat(pvalDat>0.05) = nan%0.05+eps;
    topoplot(pvalDat,LIMO.data.chanlocs,'style',style,'maplimits',[-3 log10(0.05)],...
        'electrodes','off','emarker2',emarker2,'numcontour',[linspace(-3,-1.301,4)],'whitebk','on',...
        'nosedir','+Y','headrad',0.8,'hcolor',hcolor,'plotrad',0.9);
    else
           topoplot(pvalDat,LIMO.data.chanlocs,'style','map','maplimits',[-3 log10(0.05)],...
        'electrodes','off','emarker2',emarker2,'whitebk','on',...
        'nosedir','+Y','headrad',0.8,'hcolor',hcolor,'plotrad',0.9);
    end
    caxis([-3 log10(0.05)+eps])
    if ~strcmp(g.style,'onlyCurves')
        
        colormap([cMap.p;1 1 1]),
        freezeColors;
    end
    
    
    if k ~= length(topoTimes)-1
        nextplot('byrow','delta',[-s.topo,1])
    elseif ~strcmp(g.style,'topoPlotOnly')
        
        c = cbar;
        pvalList = unique(sort([0.05 0.01 0.001 0.0001 ]));
        yTickManual = sort(log10(pvalList));
        set(c,'YTick',yTickManual)
        set(c,'YTickLabel',pvalList)
        
        newPos = get(c,'Position');
        newPos(1) = newPos(1)-0.018;
        set(c,'Position',newPos)
        xlabel(c,'min. P-Value')
        freezeColors(c)
        
    end
end
% colormap(cbrewer('div','RdYlBu',size(jet,1),[],1))
set(gcf,'Color','w')
