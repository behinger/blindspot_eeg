function bs_limo_display_constant(PathName,sigElec,varargin)

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
pred = size(LI.design.X,2);
predOrg = pred;

    g.title = 'Constant';
    

%% Define Predictors
if ischar(g)
    error(g)
end
choice = 'use theoretical p values';
LI.design.XDescAll = bs_limo_designMat_betterDesc(LI.design.XDesc);
predName = 'constant'

%% Load the Data


% We need to check whether some of the predictors are actually coded
% double, then we identify them and just sum them as everything in a GLM is
% linear anyway.
rawDataS = bs_limoGatherData('paths',LIMO.data.fullPath,'predictor',pred);

timevect = linspace(LIMO.data.start*1000,LIMO.data.end*1000,400);

%% Load the statistics

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

    dataToPlot = winMean(rawDataS(:,:,:),3);



if ~strcmp(g.style,'topoPlotOnly')
    
    %% Plot the ERP at the siginficant electrods
    
    nextplot('newRow','size',[14 s.line])
    
    h = plot(timevect,dataToPlot','Color',[0.6 0.6 0.6]);
    
    if ~isempty(g.selChan)
        set(h(g.selChan),'Color',[1 0 0])
    end
    hold all;
    
    %% Mark all alpha < 0.05
    
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
    
end


%% Plot the Topoplots

topoDiff = squeeze(mean(rawDataS(:,:,:),3));
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
    
    sigElecs = [];
    
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
        'electrodes','off','emarker2',emarker2,...
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
%     pvalDat = min(M(:,fromT:toT),[],2);
%     pvalDat = log10(pvalDat);
%     pvalDat(pvalDat==-Inf) = -4;
    %     pvalDat(pvalDat>0.05) = nan%0.05+eps;
%     topoplot(pvalDat,LIMO.data.chanlocs,'style',style,'maplimits',[-3 log10(0.05)],...
%         'electrodes','off','emarker2',emarker2,...
%         'nosedir','+Y','headrad',0.8,'hcolor',hcolor,'plotrad',0.9);
%     caxis([-3 log10(0.05)+eps])
%     if ~strcmp(g.style,'onlyCurves')
        
%         colormap([cMap.p;1 1 1]),
%         freezeColors;
%     end
    
    
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
