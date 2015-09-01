function bs_limo_display_results(Type,FileName,PathName,p,MCC,LIMO,varargin)

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


cd(PathName)
load (FileName);
% LIMO.design.bootstrap = 600;
if nargin <= 6
    flag = 1;
end
g = be_inputcheck(varargin,...
    {'interactive','string',{'yes','no'},'no';
    'title','string','','';
    'addTopo','string',{'yes','no'},'no';
    'caxis','integer',[],[];
    'elec','integer',[],[];
    'recalcMCC','boolean',[],0;
    'contour','string',{'yes','no'},'yes';
    'elecResort','string',{'yes','no'},'no';
    'topoTime','integer',[],[]});

if ischar(g)
    error(g)
end
choice = 'use theoretical p values'; % threshold based on what is computed since H0 is used for clustering
% see limo_stat_values

if LIMO.design.bootstrap == 0
    if MCC == 2 || MCC == 3
        errordlg2('Clustering thresholding necessitates boostrap - invalid choice');
    elseif MCC == 4
        errordlg2('Maximum stat thresholding necessitates bootstrap - invalid choice');
    elseif MCC == 5
        errordlg2('TFCE thresholding necessitates boostrap - invalid choice');
    end
    MCC = 1;
end



if ~strncmp(FileName,'LIMO',4) % in all cases but ERP plot for GLM
    if exist(sprintf('H0/mcc_%i_%s',MCC,FileName),'file') && ~g.recalcMCC
        load(sprintf('H0/mcc_%i_%s',MCC,FileName))
    end
    if ~exist('mask','var') || isempty(mask)
    
        [M, mask, mytitle] = limo_stat_values(Type,FileName,p,MCC,LIMO,choice,[]);
        if isnan(M)
            M = nan(size(mask));
        end
        
        
        if MCC>1 && ~isempty(mask)
            save(sprintf('H0/mcc_%i_%s',MCC,FileName),'mask','M','mytitle')
        end
        
    end
%     if isempty(mask)
%         return
%     elseif sum(mask(:)) == 0
%         warndlg('  no values under threshold  ','no significant effect','modal');
%         toplot = []; return
%     else
%         assignin('base','p_values',M)
%         assignin('base','mask',mask)
%     end
    
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
end

% ------------------------------
%      Image and topoplot
% ----------------------------
if Type == 1 || Type == 2
    
    if Type == 1
        %--------------------------
        % imagesc of the results
        %--------------------------
        % what to plot
        scale = toplot.*mask;
        v = max(scale(:));
        [e,f]=find(scale==v);
        
        %             if min(scale(:))<0
        %                 scale(scale==0)=min(scale(:))+(min(scale(:))/10);
        %             else
        %                 scale(scale==0)=NaN;
        %             end
        %             if ~isempty(g.caxis)
        %              scale = g.caxis;
        %             end
        % do the figure
        figure; set(gcf,'Color','w');
        tmp = get(gcf,'Position');
        set(gcf,'Position',[tmp(1:2) 900 500])
        ax(1) = subplot(3,3,[1 2 4 5 7 8]);
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
        
        
        
        if MCC ~=1
            [posclusterslabelmat,nposclusters] = limo_ft_findcluster(pMatrix<=p,LIMO.data.neighbouring_matrix,2);
            posclusterslabelmat = posclusterslabelmat(elecSort,:);
            alphaMask = mask(elecSort,:);
            toplotSort = toplot(elecSort,:);
            %             alphaMask(posclusterslabelmat>0) = 0.1;
            alphaMask(mask(elecSort,:)) = 1;
            
            % calculate better
            
            
            imagesc(timevect,[1:64],toplotSort,'AlphaData',alphaMask,'HitTest', 'off');
            cMap = cbrewer('div','RdYlBu',size(jet,1));
            colormap(gca,cMap(end:-1:1,:))
            hold on,
            cMap = cbrewer('qual','Accent',nposclusters);
            legendEntry = [];
            hContour = [];
            for pCl = 1:nposclusters
                
                
                %                 annotation
                [r c] = find(posclusterslabelmat==pCl,1,'last');
                p_vals(pCl) = M(elecSort(r),c);
                
                if p_vals(pCl)> 0.5 || isnan(p_vals(pCl));
                    continue
                end
                [~,tmp] = contour(timevect,[1:64],posclusterslabelmat==pCl,1, 'Edgecolor',cMap(pCl,:),'LineWidth',2);
                hContour = [hContour tmp];
                set(gca,'YDir','reverse')
                posAx = get(gca,'Position');
                %                 set(gca,'Units','pixels')
                %                 calcPos = [(r./size(toplot,1))*posAx(3) (c./size(toplot,2))*posAx(4)];
                %                 calcPos = [r c];
                %                 calcPos = [(r./size(toplot,1)) (c./size(toplot,2)) 0.1 0.1];
                %                 text(timevect(c),r,sprintf('%.4f',p_vals(pCl)), 'Color',cMap(pCl,:))
                legendEntry = [legendEntry {sprintf('%.4f',p_vals(pCl))}];
            end
            if ~isempty(hContour);      legend(hContour,legendEntry,'Location','SouthWest');end
        else
            imagesc(timevect,1:64,toplot,'HitTest', 'off');
            cMap = cbrewer('div','RdYlBu',size(jet,1));
            colormap(gca,cMap(end:-1:1,:))
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
        
        color_images_(scale,LIMO,strcmp(g.elecResort,'yes'));
        set(gca,'layer','top');
        
        if ~isempty(g.caxis)
            caxis(g.caxis) % XXX Behinger
        end
        
        if size(toplot,1)>1
            ax(2) = axes('Position',[0.56 0.1 0.5 0.5]);
            if length(f)>0
                f=f(1);
            end
            if ~isempty(g.topoTime)
                [~,f] = min(abs(timevect-g.topoTime));
            else %if isempty
                v = max(toplot(:)); [~,f]=find(toplot==v);

            end
            topoplot(toplot(:,f),LIMO.data.chanlocs);
            if ~isempty(g.caxis)
                caxis(g.caxis) % XXX Behinger
            end
            g.caxis = caxis; % save it to have the same colorsacle always!
            hCbar = cbar;
            tmp = get(hCbar,'Position');
            set(hCbar,'Position',[tmp(1) - 0.13 tmp(2:4)])
            axes(ax(2))
            %                 title(['topoplot @' num2str(timevect(f)) 'ms'],'FontSize',12)
            axTitle = title(sprintf('tvalue @ %ims',timevect(f)))
            if strcmp(g.addTopo,'yes')
                [addA] = plot_add_topo(f,LIMO,FileName);
            end
        end
        
        
        if size(toplot,1)>1
            %             update = 0;
            %             while update ==0 && strcmp(g.interactive,'yes')
            if ~exist('addA','var')
                addA = [];
            end
            %                 drawnow
            set(ax(1),'ButtonDownFcn',{@update_topo,ax(1),ax(2),LIMO,g,axTitle,addA,FileName})
            %                 set(ax(1),'Hit
            %
            %                 axes
            %
            %                 try % use try so that if figure deleted no complain
            %                     [x,y,button]=ginput(1);
            %                 catch
            %                     update = 1; break
            %                 end
            %                 if button > 1
            %                     update = 1;
            %                 end
            %                 clickedAx = gca;
            %                 if clickedAx ~=ax(1)
            %                     disp('right click to exit')
            %                 else
            %                     frame = frame_zeros + round(x / ratio);
            %                     cla(ax(2))
            %                     frame(frame < 0) = 0;
            %
            %
            %                     ax(2) =  axes('Position',[0.56 0.1 0.5 0.5]);
            %                     topoplot(toplot(:,frame),LIMO.data.chanlocs);
            %
            %                     caxis(g.caxis) % XXX Behinger
            %
            %                     drawnow
            %                     %                     cla(axTitle);
            %                     set(axTitle,'String',sprintf('tvalue @ %ims',round(x)))
            %                     %                     axTitle = title();
            %                     clear x y button
            %                     %                     tmpAx = axes;
            %                     %                     axes(axes)
            %                     %                     cla(gca)
            %                     if strcmp(g.addTopo,'yes')
            %                         [addA] = plot_add_topo(frame,LIMO,FileName,addA);
            %                     end
            %                 end
            %             end
        end
        
    elseif Type == 2
        %--------------------------
        % topoplot
        %--------------------------
        if ~isempty(LIMO.design.electrode)  % not full scalp
            msgbox('Only one electrode found','No topoplot')
        elseif sum(mask(:)) == 0
            warndlg('no values under threshold','no significant effect');
        else
            EEG.data = toplot;
            EEG.setname = mytitle;
            EEG.pnts = size(EEG.data,2);
            EEG.xmin = LIMO.data.start;
            EEG.xmax = LIMO.data.end;
            EEG.times =  (LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000); % in sec;
            EEG.trials = 1;
            EEG.chanlocs = LIMO.data.chanlocs;
            EEG.nbchan = size(EEG.data,1);
            pop_topoplot(EEG);
            assignin('base','Plotted_data',EEG.data)
        end
    end
    
    
elseif Type == 3
    
    %--------------------------
    % ERP
    %--------------------------
    
    if strncmp(FileName,'one_sample',10) || strncmp(FileName,'two_samples',11) || strncmp(FileName,'paired_samples',14) || ...
            strncmp(FileName,'con_',4) || strncmp(FileName,'ess_',4);
        % ------------------------------------------------------------------------------------------------------------
        % stat file dim = (electrodes, frames, [mean value, se, df, t, p])
        % H0 file dim = (electrodes,frames,[t, p],nboot)
        
        if strncmp(FileName,'one_sample',10); [e,f,d]=size(one_sample); data = one_sample;
        elseif strncmp(FileName,'two_samples',11); [e,f,d]=size(two_samples); data = two_samples;
        elseif  strncmp(FileName,'paired_samples',14);  [e,f,d]=size(paired_samples); data = paired_samples;
        elseif  strncmp(FileName,'con',4);  [e,f,d]=size(con); data = con;
        elseif  strncmp(FileName,'ess',4);  [e,f,d]=size(ess); data = ess;
        end
        
        if e > 1
            electrode = inputdlg('which electrode to plot','Plotting option');
            if isempty(electrode) || strcmp(cell2mat(electrode),'')
                [v,e] = max(data(:,:,4)); [v,c]=max(v); electrode = e(c);
            else
                electrode = eval(cell2mat(electrode));
                if length(electrode) > 1
                    error('1 electrode only can be plotted')
                elseif electrode > e
                    error('electrode number invalid')
                end
            end
        else
            electrode = 1;
        end
        
        timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
        trimci(:,:,2) = squeeze(data(electrode,:,1));
        if strncmp(FileName,'ess',4)
            start_at = max(strfind(FileName,'_'))+1;
            C = LIMO.contrast{eval(FileName(start_at:end-4))}.C;
            df = rank(C); % rank of the relevant contrast
            trimci(:,:,1) = squeeze(trimci(:,:,2))-finv(p./2*size(C,1),df,squeeze(data(electrode,:,3))).*squeeze(data(electrode,:,2));
            trimci(:,:,3) = squeeze(trimci(:,:,2))+finv(p./2*size(C,1),df,squeeze(data(electrode,:,3))).*squeeze(data(electrode,:,2));
        else
            trimci(:,:,1) = squeeze(trimci(:,:,2))-tinv(p./2,squeeze(data(electrode,:,3))).*squeeze(data(electrode,:,2));
            trimci(:,:,3) = squeeze(trimci(:,:,2))+tinv(p./2,squeeze(data(electrode,:,3))).*squeeze(data(electrode,:,2));
        end
        figure;set(gcf,'Color','w')
        plot(timevect,squeeze(trimci(:,:,2)),'LineWidth',3);
        fillhandle = patch([timevect fliplr(timevect)], [squeeze(trimci(:,:,1)),fliplr(squeeze(trimci(:,:,3)))], [1 0 0]);
        set(fillhandle,'EdgeColor',[1 0 1],'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
        grid on; box on; axis tight
        sig = single(mask(electrode,:)); sig(find(sig==0)) = NaN;
        h = axis;  hold on; plot(timevect,sig.*h(3),'r.','MarkerSize',20)
        set(gca,'FontSize',14,'layer','top');
        xlabel('Time in ms','FontSize',14)
        ylabel('Amplitude in {\mu}V per std','FontSize',14)
        if e>1
            title(sprintf('%s \n electrode %s (%g)',mytitle,LIMO.data.chanlocs(electrode).labels,electrode),'FontSize',16); drawnow;
        else
            title(sprintf('%s \n optimized electrode',mytitle),'FontSize',16); drawnow;
        end
        assignin('base','Plotted_data',trimci);
        
        
    elseif strncmp(LIMO.design.name,'regression analysis',19) || strncmp(LIMO.design.name,'ANOVA',5) || strncmp(LIMO.design.name,'ANCOVA',6)
        % --------------------------------------------------------------------------------------------------------------------------------
        
        % which variable(s) to plot
        % ----------------------
        if size(LIMO.design.X,2) > 2
            input_title = sprintf('which regressor to plot?: 1 to %g ',size(LIMO.design.X,2));
            regressor = inputdlg(input_title,'Plotting option'); if isempty(regressor); return; end
            try regressor = sort(eval(cell2mat(regressor)));
                if max(regressor) > size(LIMO.design.X,2); errordlg('invalid regressor number'); end
            catch ME
                return
            end
        else
            regressor = 1;
        end
        
        categorical = sum(LIMO.design.nb_conditions) + sum(LIMO.design.nb_interactions);
        if max(regressor) == size(LIMO.design.X,2); tmp = regressor(1:end-1); else tmp = regressor; end
        cat = sum(tmp<=categorical); cont = sum(tmp>categorical);
        if cat >=1 && cont >=1
            errordlg('you can''t plot categorical and continuous regressors together'); return
        end
        
        % which ERP to make
        % ------------------
        extra = questdlg('Plotting ERP','ERP Options','Original','Modelled','Adjusted','Adjusted');
        if isempty(extra)
            return;
        elseif strcmp(extra,'Original')
            if regressor == size(LIMO.design.X,2)
                errordlg('you can''t plot adjusted mean for original data'); return
            end
        end
        
        % which electrode to plot
        % ------------------
        if isempty(LIMO.design.electrode)
            electrode = inputdlg('which electrode to plot','Plotting option');
            try
                electrode = eval(cell2mat(electrode));
                if size(electrode) > 1
                    errordlg('invalid electrode number'); return
                elseif electrode > size(LIMO.data.chanlocs,2)
                    errordlg('invalid electrode number'); return
                end
            catch ME
                load R2; [v,e] = max(R2(:,:,1)); [v,c]=max(v); electrode = e(c);
                % [v,electrode] = max(max(mean(Yr,3),[],2)); % plot at the max of Yr
            end
        else
            electrode =1;  % for optimized electrode analyses
        end
        % timing info
        timevect = LIMO.data.start*1000:(1000/LIMO.data.sampling_rate):LIMO.data.end*1000; % in sec
        
        % down to business
        % ----------------------
        probs = [p/2; 1-p/2];
        z = norminv(probs);
        
        if strcmp(extra,'Original')
            load Yr;
            if sum(regressor <= categorical) == length(regressor) % for categorical variables
                for i=1:length(regressor)
                    index{i} = find(LIMO.design.X(:,regressor(i)));
                    data = squeeze(Yr(electrode,:,index{i}));
                    average(i,:) = nanmean(data,2);
                    se = (nanstd(data') ./ sqrt(numel(index{i})));
                    ci(i,:,:) = repmat(average(i,:),2,1) + repmat(se,2,1).*repmat(z,1,size(Yr,2));
                end
                if isempty(LIMO.design.electrode)
                    mytitle = sprintf('Original subjects'' parameters at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                else
                    mytitle = sprintf('Original subjects'' parameters at optimized electrode');
                end
            else % continuous variable
                for i=1:length(regressor)
                    index{i} = find(LIMO.design.X(:,regressor(i)));
                    [reg_values(i,:),sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                    continuous(i,:,:) = Yr(electrode,:,sorting_values);
                    if isempty(LIMO.design.electrode)
                        mytitle{i} = sprintf('Original subjects'' parameters \n sorted by regressor %g electrode %s (%g)', regressor(i), LIMO.data.chanlocs(electrode).labels, electrode);
                    else
                        mytitle{i} = sprintf('Original subjects'' parameters, \n sorted by regressor %g at optimized electrode', regressor(i));
                    end
                end
            end
            clear Yr
        elseif strcmp(extra,'Modelled')
            load Betas; Betas = squeeze(Betas(electrode,:,:));
            Yh = (LIMO.design.X*Betas')'; % modelled data
            load Yr; R = eye(size(Yr,3)) - (LIMO.design.X*pinv(LIMO.design.X));
            
            if sum(regressor <= categorical) == length(regressor) % for categorical variables
                for i=1:length(regressor)
                    index{i} = find(LIMO.design.X(:,regressor(i)));
                    data = squeeze(Yh(:,index{i}));
                    average(i,:) = nanmean(data,2);
                    index{i} = index{i}(find(~isnan(squeeze(Yr(electrode,1,index{i})))));
                    var   = diag(((R(index{i},index{i})*squeeze(Yr(electrode,:,index{i}))')'*(R(index{i},index{i})*squeeze(Yr(electrode,:,index{i}))')) / LIMO.model.model_df(2));
                    CI = sqrt(var/size(index{i},1))*z';
                    ci(i,:,:) = (repmat(nanmean(data,2),1,2)+CI)';
                end
                if isempty(LIMO.design.electrode)
                    mytitle = sprintf('Modelled subjects'' parameters at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                else
                    mytitle = sprintf('Modelled subjects'' parameters at optimized electrode');
                end
            else % continuous variable
                for i=1:length(regressor)
                    index{i} = find(LIMO.design.X(:,regressor(i)));
                    [reg_values(i,:),sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                    continuous(i,:,:) = Yh(:,sorting_values);
                    if isempty(LIMO.design.electrode)
                        mytitle{i} = sprintf('Modelled single subjects'' parameters \n sorted by regressor %g electrode %s (%g)', regressor(i), LIMO.data.chanlocs(electrode).labels, electrode);
                    else
                        mytitle{i} = sprintf('Modelled single subjects'' parameters \n sorted by regressor %g at optimized electrode', regressor(i));
                    end
                end
            end
        else % Adjusted
            all = [1:size(LIMO.design.X,2)-1]; all(regressor)=[];
            load Yr; Yr = squeeze(Yr(electrode,:,:));
            load Betas; Betas = squeeze(Betas(electrode,:,:));
            confounds = (LIMO.design.X(:,all)*Betas(:,all)')';
            Ya = Yr - confounds; clear Yr Betas confounds;
            if sum(regressor <= categorical) == length(regressor) % for categorical variables
                for i=1:length(regressor)
                    index{i} = find(LIMO.design.X(:,regressor(i)));
                    data = squeeze(Ya(:,index{i}));
                    average(i,:) = nanmean(data,2);
                    se = nanstd(data') ./ sqrt(numel(index{i}));
                    ci(i,:,:) = repmat(average(i,:),2,1) + repmat(se,2,1).*repmat(z,1,size(Ya,1));
                end
                if isempty(LIMO.design.electrode)
                    mytitle = sprintf('Adjusted subjects'' parameters at electrode %s (%g)', LIMO.data.chanlocs(electrode).labels, electrode);
                else
                    mytitle = sprintf('Adjusted subjects'' parameters at  at optimized electrode');
                end
            else % continuous variable
                for i=1:length(regressor)
                    index{i} = find(LIMO.design.X(:,regressor(i)));
                    [reg_values(i,:),sorting_values]=sort(LIMO.design.X(index{i},regressor(i)));  % continuous variable 3D plot
                    continuous(i,:,:) = Ya(:,sorting_values);
                    if isempty(LIMO.design.electrode)
                        mytitle{i} = sprintf('Adjusted single subjects'' parameters \n sorted by regressor %g electrode %s (%g)', regressor(i), LIMO.data.chanlocs(electrode).labels, electrode);
                    else
                        mytitle{i} = sprintf('Adjusted single subjects'' parameters \n sorted by regressor %g at optimized electrode', regressor(i));
                    end
                end
            end
        end
        
        % make the figure(s)
        % ------------------
        figure;set(gcf,'Color','w')
        if sum(regressor <= categorical) == length(regressor)
            for i=1:size(average,1)
                plot(timevect,squeeze(average(i,:)),'LineWidth',1.5); hold on
                if i==1
                    colorOrder = get(gca, 'ColorOrder');
                    colorOrder = repmat(colorOrder,ceil(size(average,1)/size(colorOrder,1)),1);
                end
                x = squeeze(ci(i,1,:)); y = squeeze(ci(i,2,:));
                fillhandle = patch([timevect fliplr(timevect)], [x',fliplr(y')], colorOrder(i,:));
                set(fillhandle,'EdgeColor',colorOrder(i,:),'FaceAlpha',0.2,'EdgeAlpha',0.8);%set edge color
            end
            
            % if regressor spans columns of an effect, plot significant time frames
            index = 1; index2 = LIMO.design.nb_conditions(1);
            for i=1:length(LIMO.design.nb_conditions)
                effect = index:index2;
                if length(regressor) == length(effect)
                    if mean(regressor == effect) == 1
                        name = sprintf('Condition_effect_%g.mat',i);
                        load(name); [M, mask, mytitle2] = limo_stat_values(1,name,p,MCC,LIMO,choice);
                        sig = single(mask(electrode,:)); sig(find(sig==0)) = NaN;
                        h = axis;  plot(timevect,sig.*h(3),'r*','LineWidth',2)
                        break
                    end
                else
                    index = index+LIMO.design.nb_conditions(i);
                    if i<length(LIMO.design.nb_conditions)
                        index2 = LIMO.design.nb_conditions(i)+LIMO.design.nb_conditions(i+1);
                    end
                end
            end
            
            if LIMO.design.nb_interactions ~= 0
                index = sum(LIMO.design.nb_conditions)+1; index2 = sum(LIMO.design.nb_conditions)+LIMO.design.nb_interactions(1);
                for i=1:length(LIMO.design.nb_interactions)
                    effect = index:index2;
                    if length(regressor) == length(effect)
                        if mean(regressor == effect) == 1
                            name = sprintf('Interaction_effect_%g.mat',i);
                            load(name); [M, mask, mytitle2] = limo_stat_values(1,name,p,MCC,LIMO,choice);
                            sig = single(mask(electrode,:)); sig(find(sig==0)) = NaN;
                            h = axis;  plot(timevect,sig.*h(3),'r*','LineWidth',2)
                            break
                        end
                    else
                        index = index+LIMO.design.nb_interactions(i);
                        if i<length(LIMO.design.nb_interactions)
                            index2 = LIMO.design.nb_interactions(i)+LIMO.design.nb_interactions(i+1);
                        end
                    end
                end
            end
            
            
            % --
            axis tight; grid on; box on
            title(mytitle,'FontSize',16); drawnow;
            assignin('base','Plotted_data', average)
            v=axis;axis([v(1) v(2) v(3)+.1*v(3) v(4)+.1*v(4)])
            set(gca,'FontSize',14);
            ylabel('Amplitude in {\mu}V','FontSize',16)
            xlabel('Time in ms','FontSize',16)
            
        else
            for i=1:size(continuous,1)
                if i > 1; figure;set(gcf,'Color','w'); end
                index = find(~isnan(squeeze(continuous(i,1,:))));
                surf(index,timevect,squeeze(continuous(i,:,index)));shading interp
                
                % --
                axis tight; title(mytitle{i},'FontSize',14); drawnow;
                xlabel('Sorted subjects','FontSize',16)
                try
                    set(gca,'XTick',index, 'XTickLabels', reg_values(index));
                end
                ylabel('Time in ms','FontSize',16)
                zlabel('Amplitude in {\mu}V per std','FontSize',16)
            end
        end
        
    else
        errordlg('this file is not supported for this kind of plot','Nothing plotted')
    end
end % closes type

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