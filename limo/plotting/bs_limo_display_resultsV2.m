function bs_limo_display_resultsV2(pred,PathName,varargin)

% This function displays various results
% The arguments specify cases for the
% different kind of figures, thresholds etc ..
%
% FORMAT:
% limo_display_results(Type,FileName,PathName,p,LIMO)
%
% INPUTS:
%   Type      = type of images/plot to do
%               1 - 2D images with a intensity plotted as function of time (x) and electrodes (y)
%               2 - topographic plot a la eeglab
%               3 - plot the ERP data (original or modeled)
%   Filename  = Name of the file to image
%   PathName  = Path of the file to image
%   p         = threshold p value e.g. 0.05
%   MCC       = Multiple Comparison technique
%               1=None, 2=2D Cluster, 3=1D Cluster, 4=T max, 5=TFCE
%   LIMO      = LIMO structure
%   flag      = indicates to allow surfing the  (1) or not (0)
%
% Cyril Pernet, Guillaume Rousselet v3 06-05-2009
% Carl Gaspar 03-09-2009 - fixed some axis issues for 3D plots (see subfunction time_vect_)
% Cyril P. v4 09-09-2009 allows random effect results to be displayed (+ some clean up)
% Cyril P. v5. 10-03-2010 split the whole file into 2 parts based on LIMO.level (1 or 2)
% Guillaume Rousselet v4 06-07-2010 added the max(T)/max(F) and cluster stats for random effect
% Cyril Pernet v4 16-05-2010 fixed the random effect to automatically load bootstrap and get the neighbouring matrix for clusters
% Nicolas Chauveau 08-12-2011 fixed the ERP plot of gp*repeated measures (for levels>2)
% Cyril Pernet v5 10-10-2012 added tfce and redesigned CI with filling
% -----------------------------
%  Copyright (C) LIMO Team 2010



% LIMO.design.bootstrap = 600;

g = be_inputcheck(varargin,...
    {'interactive','string',{'yes','no'},'no';
    'title','string','','';
    'addTopo','string',{'yes','no'},'no';
    'caxis','integer',[],[];
    'elec','integer',[],[];
    'recalcMCC','boolean',[],0;
    'contour','string',{'yes','no'},'yes';
    'elecResort','string',{'yes','no'},'no';
    'MCC','integer',[1:5],2;
    'alpha','real',[0 1],0.05;
    'interaction','string',{'yes','no'},'no';
    'topoTime','integer',[],[]});


cd(PathName)
load('LIMO.mat')
[~,LI] = bs_limo_plotDesignmatrix(LIMO);
[~,a] = fileparts(LIMO.dir);
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

if ischar(g)
    error(g)
end
choice = 'use theoretical p values';
%% Load The clustering mask
%% We want to show Factor-Values, T-values (later Semi-Partial
%% Correlations) and Topo-Difference activations.
if strncmp(FileName,'R2',2)
    toplot = squeeze(R2(:,:,1));
    assignin('base','R_values',toplot)
elseif strncmp(FileName,'one_sample',10)
    toplot = squeeze(one_sample(:,:,4));
    pMatrix =  squeeze(one_sample(:,:,5));
    
    assignin('base','T_values',toplot)
elseif strncmp(FileName,'two_samples',11)
    toplot = squeeze(two_samples(:,:,4));
    assignin('base','T_values',toplot)
    
    
elseif strncmp(FileName,'paired',6)
    toplot = squeeze(paired_samples(:,:,4));
    pMatrix =  squeeze(paired_samples(:,:,5));
    assignin('base','T_values',toplot)
    
    
    
elseif strncmp(FileName,'Covariate',9);
    toplot = squeeze(Covariate_effect(:,:,1));
    assignin('base','F_values',toplot)
    
elseif strncmp(FileName,'Condition',9);
    toplot = sqrt(squeeze(Condition_effect(:,:,1)));
    pMatrix =  squeeze(Condition_effect(:,:,2));
    assignin('base','F_values',toplot)
    fprintf('only F-Values, no t-values available, showing sqrt(F-Value)')
    
elseif strncmp(FileName,'Interaction',9);
    toplot = sqrt(squeeze(Interaction_effect(:,:,1)));
    pMatrix =  squeeze(Interaction_effect(:,:,2));
    assignin('base','F_values',toplot)
    g.addTopo = 'no';
    fprintf('deactivated topoplots cause of interactions')
    
elseif strncmp(FileName,'con_',4);
    toplot = squeeze(con(:,:,4));
    assignin('base','T_values',toplot)
elseif strncmp(FileName,'ess_',4);
    toplot = squeeze(ess(:,:,4));
    assignin('base','F_values',toplot)
elseif strncmp(FileName,'Rep_ANOVA_Interaction',21);
    toplot = squeeze(Rep_ANOVA_Interaction_with_gp(:,:,1));
    assignin('base','F_values',toplot)
elseif strncmp(FileName,'Rep_ANOVA_Gp',12);
    toplot = squeeze(Rep_ANOVA_Gp_effect(:,:,1));
    assignin('base','F_values',toplot)
elseif strncmp(FileName,'Rep_ANOVA',9);
    toplot = squeeze(Rep_ANOVA(:,:,1));
    assignin('base','F_values',toplot)
    pMatrix =  squeeze(Rep_ANOVA(:,:,2));
    
end

%%
if g.MCC == 5
[mask M] = bs_limo_tfceCalc(g.MCC,FileName,g.recalcMCC,LIMO,toplot);
elseif g.MCC == 2
  if exist(sprintf('H0/mcc_%i_%s',g.MCC,FileName),'file') && ~g.recalcMCC
    load(sprintf('H0/mcc_%i_%s',g.MCC,FileName))
   end 
end


if ~exist('mask','var') || isempty(mask)
    Type = 1
    [M, mask, mytitle] = limo_stat_values(Type,FileName,g.alpha,g.MCC,LIMO,choice,[]);
    if isnan(M)
        M = nan(size(mask));
    end
    
    
    if g.MCC>1 && ~isempty(mask)
        save(sprintf('H0/mcc_%i_%s',g.MCC,FileName),'mask','M','mytitle')
    end
    
end




% ------------------------------
%      Image and topoplot
% ----------------------------


%--------------------------
% imagesc of the results
%--------------------------
% what to plot
scale = toplot.*mask;
v = max(scale(:));
[e,f]=find(scale==v);


% do the figure
figure; set(gcf,'Color','w');
tmp = get(gcf,'Position');
set(gcf,'Position',[tmp(1:2) 1400 800])
% ax(1) = subplot(3,3,[1 2 4 5 7 8]);

jisubplot(20,20)
nextplot('size',[12,9])
timevect = linspace(LIMO.data.start*1000,LIMO.data.end*1000,size(toplot,2));

ratio = (LIMO.data.end*1000 - LIMO.data.start*1000) / size(toplot,2);
if LIMO.data.start < 0
    frame_zeros = round(abs(LIMO.data.start*1000) / ratio);
end

if strcmp(g.elecResort,'yes')
    elecSort = bs_elecResort(LIMO.data.chanlocs);
    LIMO.data.elecSort = elecSort;
else
    elecSort = 1:size(toplot,1);
end


toplotSort = toplot(elecSort,:);

if g.MCC  == 2
    [posclusterslabelmat,nposclusters] = limo_ft_findcluster(pMatrix<=g.alpha,LIMO.data.neighbouring_matrix,2);
    posclusterslabelmat = posclusterslabelmat(elecSort,:);
    alphaMask = mask(elecSort,:);
    %             alphaMask(posclusterslabelmat>0) = 0.1;
    alphaMask(mask(elecSort,:)) = 1;
    
    imagesc(timevect,[1:64],toplotSort,'AlphaData',alphaMask,'HitTest', 'off');
    cMap = cbrewer('div','RdYlBu',size(jet,1));
    colormap(gca,cMap(end:-1:1,:))
    hold on,
    cMap = cbrewer('qual','Accent',nposclusters);
    legendEntry = [];
    hContour = [];
    for pCl = 1:nposclusters
        [r c] = find(posclusterslabelmat==pCl,1,'last');
        p_vals(pCl) = M(elecSort(r),c);
        if p_vals(pCl)> 0.5 || isnan(p_vals(pCl));
            continue
        end
        [~,tmp] = contour(timevect,[1:64],posclusterslabelmat==pCl,1, 'Edgecolor',cMap(pCl,:),'LineWidth',2);
        hContour = [hContour tmp];
        set(gca,'YDir','reverse')
        posAx = get(gca,'Position');
        legendEntry = [legendEntry {sprintf('%.4f',p_vals(pCl))}];
    end
    if ~isempty(hContour);      legend(hContour,legendEntry,'Location','SouthWest');end
    color_images_(scale,LIMO,strcmp(g.elecResort,'yes'));
elseif g.MCC == 5
    alphaMask = mask(elecSort,:);
    Mresort = M(elecSort,:);
    imagesc(timevect,[1:64],log10(Mresort),'AlphaData',alphaMask,'HitTest', 'off');
    colormap([cbrewer('seq','RdPu',100,[],1)]);
    caxis([-3 log10(0.05)+eps])
    
    freezeColors
    
    h = colorbar('location','West');
    y = get(h,'YTick');
    set(h,'YTickLabel',round((10.^y)*1000)/1000);
    cbfreeze
    
    
else
    imagesc(timevect,1:64,toplotSort,'HitTest', 'off');
    cMap = cbrewer('div','RdYlBu',size(jet,1));
    color_images_(scale,LIMO,strcmp(g.elecResort,'yes'));
end
%             figure,imagesc(timevect,1:size(toplot,1),posclusterslabelmat);
%             aContour = axes
%             set(aContour,           'Color','none','Position',get(ax(1),'Position'))
%             imagesc(timevect,1:size(toplot,1),posclusterslabelmat,'AlphaD
%             ata',alphaMask2*0.1);
%             hold on
%             contour(timevect,1:size(toplot,1),mask,[1],'k')
%             imagesc(timevect,1:size(toplot,1),mask)

title([g.title],'FontSize',16);


set(gca,'layer','top');

if ~isempty(g.caxis)
    caxis(g.caxis) % XXX Behinger
end
%%

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

rawDataS = bs_limoGatherData('paths',LIMO.data.fullPath,'predictor',preds);

if  strcmp(g.interaction,'yes')
    rawData = rawDataS(:,:,1,:) - rawDataS(:,:,2,:)- rawDataS(:,:,3,:) + rawDataS(:,:,4,:);
    
else
    rawData = rawDataS(:,:,1,:) - rawDataS(:,:,2,:);
    
end


%% Plot the ERP at the siginficant electrods
sigElecs = any(mask(:,:)');
% sigElecs = 46;
% sigElecs = 49;
sigElecs = 1:64;
nextplot('newRow')%('size',[10 10])
% tmpDat = permute(rawDataS(sigElecs,:,:,:),[2 4 1 3]);
% topoSTD = squeeze(std(rawData,[],1));

topoSTD = squeeze(diff(prctile(rawData,[5 95],1)));
% topoSTD = squeeze(max(rawData,[],1) - min(rawData,[],1));
% topoSTD = mean(topoSTD,2);

h = plot(repmat(timevect,size(topoSTD,2),1)',topoSTD);hold all;
% if strcmp(g.interaction,'yes')
%     hk = plot(repmat(timevect,4,1)',squeeze(mean(std(rawDataS(sigElecs,:,:,:),[],1),4)));
% end
hM = plot(timevect,squeeze(mean(topoSTD,2))','LineWidth',3);
[~, LI] = bs_limo_plotDesignmatrix(LIMO);

if strcmp(g.interaction,'yes')
    startIdx = length(LI.design.nb_conditions);
%     legend([hM;hk],[{'Variance over Electrodes Interaction Betas'},LI.design.XDesc(startIdx+(l-1)*4+1:startIdx+4+(l-1)*4)],'location','NorthWest');
    legend([hM],[LI.design.XDesc(startIdx+(l-1)*4+1)],'location','NorthWest');
    set(gca,'YLim',[0 10]);
else
    legend([hM],{'Variance over Electrodes of Beta'})    ;
    set(gca,'YLim',[0 10]);
end
topoDiff = squeeze(mean(rawData(sigElecs,:,:,:),1));
title('Pred1 (off) - Pred2 (On)')

interactLinePlot([hM;h])
set(gca,'Box','off')
hline(0,'--k')


%% Plot the Topoplots
% nextplot('size',[10 20])
topoDiff = squeeze(mean(rawData(:,:,:,1),4));
scale = max([abs(min(topoDiff(:))) abs(max(topoDiff(:)))]);
topoTimes = [-300:25:400];
for k = 1:length(topoTimes)-1
    if k == 1
        nextplot('size',2)
        nextplot('delta',[-8 0])
        
    elseif mod(k,4) == 1
        nextplot('delta',[2 -6])
    else
        nextplot('size',2)
        
    end
    fromT = find(timevect<=topoTimes(k),1,'last');
    toT = find(timevect<topoTimes(k+1),1,'last');
    tmpDat= mean(topoDiff(:,fromT:toT),2);
    sigElecs = any(mask(:,fromT:toT)');
    
    topoplot(tmpDat,LIMO.data.chanlocs,'maplimits',[-scale scale],'electrodes','off','emarker2',{[find(sigElecs)] '.','k',10,1},'nosedir','+Y');
    hTitle = title(sprintf('%.0f <= %.0f',timevect(fromT),timevect(toT)));
    tmpPos = get(hTitle,'Position'); tmpPos(2) = tmpPos(2)-0.3;
    set(hTitle,'Position',tmpPos);
    
    
end
colormap(cbrewer('div','RdYlBu',size(jet,1),[],1))
cbar
%%


%
% if size(toplot,1)>1
%     ax(2) = axes('Position',[0.56 0.1 0.5 0.5]);
%     if length(f)>0
%         f=f(1);
%     end
%     if ~isempty(g.topoTime)
%         [~,f] = min(abs(timevect-g.topoTime));
%     else %if isempty
%         v = max(toplot(:)); [~,f]=find(toplot==v);
%
%     end
%     topoplot(toplot(:,f),LIMO.data.chanlocs);
%     if ~isempty(g.caxis)
%         caxis(g.caxis) % XXX Behinger
%     end
%     g.caxis = caxis; % save it to have the same colorsacle always!
%     hCbar = cbar;
%     tmp = get(hCbar,'Position');
%     set(hCbar,'Position',[tmp(1) - 0.13 tmp(2:4)])
%     axes(ax(2))
%     %                 title(['topoplot @' num2str(timevect(f)) 'ms'],'FontSize',12)
%     axTitle = title(sprintf('tvalue @ %ims',timevect(f)))
%     if strcmp(g.addTopo,'yes')
%         [addA] = plot_add_topo(f,LIMO,FileName);
%     end
% end
%


end % closes the function

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                                   ROUTINES
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% color map
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function color_images_(scale,LIMO,sortEl)

% cc=colormap(jet);
% cc = colormap(ametrine);
% cc(1,:)=[.9 .9 .9];
% colormap(gca,cc);
set(gca,'XMinorTick','on','LineWidth',2)
try
    set(gca,'YTick',1:length(LIMO.data.expected_chanlocs));
catch ME
    set(gca,'YTick',1:length(LIMO.data.chanlocs));
end

ylabel('Electrodes','FontSize',14);
xlabel('Time in ms','FontSize',14);

if LIMO.Level == 1
    for i = 1 : length(LIMO.data.chanlocs)
        %         try
        if sortEl
            label_electrodes{i} = LIMO.data.chanlocs(LIMO.data.elecSort(i)).labels; %behinger
        else
            label_electrodes{i} = LIMO.data.chanlocs(i).labels; %behinger
        end
        %             label_electrodes{i} = LIMO.data.expected_chanlocs(i).labels;
        %         catch ME
        %             label_electrodes{i} = LIMO.data.chanlocs(i).labels;
        %         end
    end
    
else
    if isempty(LIMO.design.electrode)
        for i = 1 : length(LIMO.data.chanlocs)
            if sortEl
                label_electrodes{i} = LIMO.data.chanlocs(LIMO.data.elecSort(i)).labels; %behinger
            else
                label_electrodes{i} = LIMO.data.chanlocs(i).labels; %behinger
            end
        end
    else
        if length(LIMO.design.electrode) == 1
            label_electrodes = LIMO.design.electrode;
        else
            label_electrodes = ' ';
            ylabel('optimized electrode','FontSize',14);
        end
    end
end
set(gca,'YTickLabel', label_electrodes);

end


%% time vector and label
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function [timevect, label] = labels_time_(LIMO,ah)
interval = 50/(1000/LIMO.data.sampling_rate); % in frame
timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
zero_column = find(timevect == 0);
if isempty(zero_column) == 1 % in case it does not encompasses 0
    zero_column = 1;
end

if  LIMO.data.start < 0
    positive_label = timevect(zero_column:interval:end);
    negative_label = timevect(zero_column:-interval:1);
    if negative_label == 0
        negative_label = LIMO.data.start*1000;
    end
    label = [fliplr(negative_label(2:end)) positive_label];
else
    label = timevect(zero_column:interval:end);
end
set(ah,'XTick',find(timevect==label(1)):interval:find(timevect==label(end)),'XTickLabel',round(label));
end

%%
function [addA] = plot_add_topo(f,LIMO,FileName,smallAx)
if nargin==4
    try
        delete(smallAx);
    catch
        % hold off
    end
end


if length(FileName) == 39
    preds = [ str2num(FileName(32:end-6)) str2num(FileName(34:end-4))];
else
    
    preds = [str2num(FileName(end-5:end-5)) str2num(FileName(end-4:end-4))];
end
if length(preds) == 1
    preds = [preds preds+1];
end
rawData = bs_limoGatherData('paths',LIMO.data.fullPath,'predictor',preds);
%first pred
% addA(1) = axes('Position',[0.65 0.70 0.15 0.2]);
addA(1) = axes('Position',[0.61 0.70 0.15 0.15]);
topoData = mean(rawData,4); %do mean over subjects
topoplot(topoData(:,f,1),LIMO.data.chanlocs);
limAddTopop = max(max(max(topoData(:,f,:))),abs(min(min(topoData(:,f,:)))));
caxis([-limAddTopop limAddTopop])  % XXX Behinger
title(sprintf('Pred %i',preds(1)),'FontSize',12)
%second pred
% addA(2) = axes('Position',[0.83 0.70 0.15 0.2]);
addA(2) = axes('Position',[0.72 0.70 0.15 0.15]);
topoplot(topoData(:,f,2),LIMO.data.chanlocs);
caxis([-limAddTopop limAddTopop])  % XXX Behinger
title(sprintf('Pred %i',preds(2)),'FontSize',12)
%difference
% addA(3) = axes('Position',[0.85 0.45 0.15 0.2]);
addA(3) = axes('Position',[0.83 0.7 0.15 0.15]);
topoplot(topoData(:,f,1)-topoData(:,f,2),LIMO.data.chanlocs);
caxis([-limAddTopop limAddTopop])  % XXX Behinger
title(sprintf('Diff %i-%i',preds(1),preds(2)),'FontSize',12)
hCbar = cbar;
% inspect(hCbar)
tmp = get(hCbar,'Position');
set(hCbar,'Position',[tmp(1) - 0.04 tmp(2:4)])
addA(4) = hCbar;
end


function update_topo(src,eventdata,imagescAx,topAx,LIMO,g,axTitle,addA,FileName)
[frameX,frameY,toplot] = getimage(imagescAx);
x = get(gca,'CurrentPoint');
x = x(1,1);
[~,frame] = min(abs(frameX-x));

frame(frame <= 0) = 1;

cla(topAx)
axes(topAx)
topoplot(toplot(:,frame),LIMO.data.chanlocs);

caxis(g.caxis) % XXX Behinger

drawnow
%                     cla(axTitle);
title(sprintf('tvalue @ %ims',round(x)))

if strcmp(g.addTopo,'yes')
    [addA] = plot_add_topo(frame,LIMO,FileName,addA);
end
end