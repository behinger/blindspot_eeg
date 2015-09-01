function [] = bs_limo_eeg_4_parfor(eIdx)
fprintf('starting electrode fit \n')
% NBOOT (updated if specified in LIMO.design)
% ------------------------------------------
nboot =  599;
% ----------

% get the LIMO.mat
files = dir;
load_limo = 0;
for i=1:size(files,1)
    if strcmp(files(i).name,'LIMO.mat')
        load('LIMO.mat');
        load_limo = 1;
    end
end

if load_limo == 0
    [file,dir_path] = uigetfile('LIMO.mat','select a LIMO.mat file');
    if file ==0
        return
    else
        cd (dir_path); load LIMO.mat;
    end
end
cd (LIMO.dir);


% ---------------- univariate analysis ------------------
% --------------------------------------------------------

if strcmp(LIMO.design.type_of_analysis,'Mass-univariate') && ~exist('Betas','file')
    
    % --------- load files created by limo_design_matrix ------------------
    load Yr; load Yhat; load Res; load R2; load Betas;
    
    
    % ------------- prepare weight matrice  -------------------------------------
    if strcmp(LIMO.design.method,'WLS') || strcmp(LIMO.design.method,'OLS')
        W = ones(size(Yr,1),size(Yr,3));
    elseif strcmp(LIMO.design.method,'IRLS')
        W = zeros(size(Yr));
    end
    
    % ------------ prepare condition/covariates -------------------
    if LIMO.design.nb_conditions ~=0
        tmp_Condition_effect = NaN(size(Yr,1),size(Yr,2),length(LIMO.design.nb_conditions),2);
    end
    
    if LIMO.design.nb_interactions ~=0
        tmp_Interaction_effect = NaN(size(Yr,1),size(Yr,2),length(LIMO.design.nb_interactions),2);
    end
    
    if LIMO.design.nb_continuous ~=0
        tmp_Covariate_effect = NaN(size(Yr,1),size(Yr,2),LIMO.design.nb_continuous,2);
    end
    
    % -------------- loop the analysis electrode per electrode
    if size(Yr,1) == 1
        array = 1;
    else
        array = find(~isnan(Yr(:,1,1))); % skip empty electrodes
    end
    
    if strcmp(LIMO.design.status,'to do')
        update = 1;
        X = LIMO.design.X;
        for e = 1:size(array,1)
            electrode = array(e); warning off;
            fprintf('analyzing electrode %g/%g \n',electrode,size(Yr,1));
            if LIMO.Level == 2
                Y = squeeze(Yr(electrode,:,:));
                index = find(~isnan(Y(1,:)));
                Y = Y(:,index);
                LIMO.design.X = X(index,:);
                model = limo_glm1(Y',LIMO); warning on;
                if isempty(index)
                    index = [1:size(Y,2)];
                end
            else % level 1 we should not have any NaNs
                index = [1:size(Yr,3)];
                model = limo_glm1(squeeze(Yr(electrode,:,:))',LIMO);
            end
            
            % update the LIMO.mat (do it only once)
            if update == 1
                LIMO.model.model_df = model.df;
                if LIMO.design.nb_conditions ~=0
                    LIMO.model.conditions_df  = model.conditions.df;
                end
                if LIMO.design.nb_interactions ~=0
                    LIMO.model.interactions_df  = model.interactions.df;
                end
                if LIMO.design.nb_continuous ~=0
                    LIMO.model.continuous_df  = model.continuous.df;
                end
                update = 0;
            end
            
            % update the files to be stored on the disk
            if  strcmp(LIMO.design.method,'IRLS')
                W(electrode,:,index) = model.W;
            else
                W(electrode,index) = model.W;
            end
            fitted_data = LIMO.design.X*model.betas;
            Yhat(electrode,:,index) = fitted_data';
            Res(electrode,:,index)  = squeeze(Yr(electrode,:,index)) - fitted_data'; clear fitted_data
            R2(electrode,:,1) = model.R2_univariate;
            R2(electrode,:,2) = model.F;
            R2(electrode,:,3) = model.p;
            Betas(electrode,:,:,1) = model.betas';
            
            if prod(LIMO.design.nb_conditions) ~=0
                if length(LIMO.design.nb_conditions) == 1
                    tmp_Condition_effect(electrode,:,1,1) = model.conditions.F;
                    tmp_Condition_effect(electrode,:,1,2) = model.conditions.p;
                else
                    for i=1:length(LIMO.design.nb_conditions)
                        tmp_Condition_effect(electrode,:,i,1) = model.conditions.F(i,:);
                        tmp_Condition_effect(electrode,:,i,2) = model.conditions.p(i,:);
                    end
                end
            end
            
            if LIMO.design.fullfactorial == 1
                if length(LIMO.design.nb_interactions) == 1
                    tmp_Interaction_effect(electrode,:,1,1) = model.interactions.F;
                    tmp_Interaction_effect(electrode,:,1,2) = model.interactions.p;
                else
                    for i=1:length(LIMO.design.nb_interactions)
                        tmp_Interaction_effect(electrode,:,i,1) = model.interactions.F(i,:);
                        tmp_Interaction_effect(electrode,:,i,2) = model.interactions.p(i,:);
                    end
                end
            end
            
            if LIMO.design.nb_continuous ~=0
                if LIMO.design.nb_continuous == 1
                    tmp_Covariate_effect(electrode,:,1,1) = model.continuous.F;
                    tmp_Covariate_effect(electrode,:,1,2) = model.continuous.p;
                else
                    for i=1:LIMO.design.nb_continuous
                        tmp_Covariate_effect(electrode,:,i,1) = model.continuous.F(i,:);
                        tmp_Covariate_effect(electrode,:,i,2) = model.continuous.p(i,:);
                    end
                end
            end
        end
        
        % save data on the disk and clean out
        LIMO.design.X       = X;
        LIMO.design.weights = W;
        LIMO.design.status = 'done';
        save LIMO LIMO; save Yhat Yhat;
        save Res Res; save Betas Betas;
        save R2 R2; clear Yhat Res Betas R2
        
        if prod(LIMO.design.nb_conditions) ~=0
            for i=1:length(LIMO.design.nb_conditions)
                name = sprintf('Condition_effect_%g',i);
                if size(tmp_Condition_effect,1) == 1
                    tmp = squeeze(tmp_Condition_effect(1,:,i,:));
                    Condition_effect = NaN(1,size(tmp_Condition_effect,2),2);
                    Condition_effect(1,:,:) = tmp;
                else
                    Condition_effect = squeeze(tmp_Condition_effect(:,:,i,:));
                end
                save(name,'Condition_effect','-v7.3')
            end
            clear Condition_effect tmp_Condition_effect
        end
        
        if LIMO.design.fullfactorial == 1
            for i=1:length(LIMO.design.nb_interactions)
                name = sprintf('Interaction_effect_%g',i);
                if size(tmp_Interaction_effect,1) == 1
                    tmp = squeeze(tmp_Interaction_effect(1,:,i,:));
                    Interaction_effect = NaN(1,size(tmp_Interaction_effect,2),2);
                    Interaction_effect(1,:,:) = tmp;
                else
                    Interaction_effect = squeeze(tmp_Interaction_effect(:,:,i,:));
                end
                save(name,'Interaction_effect','-v7.3')
            end
            clear Interaction_effect tmp_Interaction_effect
        end
        
        if LIMO.design.nb_continuous ~=0
            for i=1:LIMO.design.nb_continuous
                name = sprintf('Covariate_effect_%g',i);
                if size(tmp_Covariate_effect,1) == 1
                    tmp = squeeze(tmp_Covariate_effect(1,:,i,:));
                    Covariate_effect = NaN(1,size(tmp_Covariate_effect,2),2);
                    Covariate_effect(1,:,:) = tmp;
                else
                    Covariate_effect = squeeze(tmp_Covariate_effect(:,:,i,:));
                end
                save(name,'Covariate_effect','-v7.3')
            end
            clear Covariate_effect tmp_Covariate_effect
        end
        clear file electrode filename model reg dir i W
    end
    fprintf('Stopped electrode fit \n Starting Bootstraping \n')
    % as above for bootstrap under H0
    % -------------------------------
%     boot_go = 0;
%     if LIMO.design.bootstrap ~=0
%         if exist('H0','dir')
%             if strcmp(questdlg('H0 directory detected, overwrite?','data check','Yes','No','No'),'No');
%                 if LIMO.design.tfce == 1
%                     errordlg2('bootstrap skipped - attempting to continue with tfce');
%                 else
%                     return
%                 end
%             else
%                 boot_go = 1;
%             end
%         else
%             boot_go = 1;
%         end
%     end
end
boot_go = 1;

    if boot_go == 1
        try
            fprintf('\n %%%%%%%%%%%%%%%%%%%%%%%% \n Bootstrapping data with the GLM can take a while, be patient .. \n %%%%%%%%%%%%%%%%%%%%%%%% \n')
            mkdir H0; load Yr;
            
            if LIMO.design.bootstrap > 599
                nboot = LIMO.design.bootstrap;
            end
            
            if ~exist('boot_table','file')
            if LIMO.Level == 2
                boot_table = limo_create_boot_table(Yr,nboot);
            else
                boot_table = randi(size(Yr,3),size(Yr,3),nboot);
            end
            save boot_table boot_table
            else
                load boot_table
            end

            H0_Betas = NaN(size(Yr,1), size(Yr,2), size(LIMO.design.X,2), nboot);
            H0_R2 = NaN(size(Yr,1), size(Yr,2), 3, nboot); % stores R, F and p values for each boot
            
            if LIMO.design.nb_conditions ~= 0
                tmp_H0_Conditions = NaN(size(Yr,1), size(Yr,2), length(LIMO.design.nb_continuous), 2, nboot);
            end
            
            if LIMO.design.nb_interactions ~=0
                tmp_H0_Interaction_effect = NaN(size(Yr,1),size(Yr,2),length(LIMO.design.nb_interactions), 2, nboot);
            end
            
            if LIMO.design.nb_continuous ~= 0
                tmp_H0_Covariates = NaN(size(Yr,1), size(Yr,2), LIMO.design.nb_continuous, 2, nboot);
            end
            
            warning off;
            W = LIMO.design.weights;
            X = LIMO.design.X;
            %                     h = waitbar(0,'bootstraping data','name','% done');
            if matlabpool('size') <= 0
                matlabpool local 7
            end
            parfor e = 1:size(array,1)
                electrode = array(e);
                %                         waitbar(e/size(array,1))
                fprintf('bootstrapping electrode %g \n',electrode);
                if LIMO.Level == 2
                    Y = squeeze(Yr(electrode,:,:));
                    index = find(~isnan(Y(1,:)));
                    if numel(size(LIMO.design.weights)) == 3
                        model = limo_glm1_boot(Y(:,index)',X(index,:),LIMO.design.nb_conditions,LIMO.design.nb_interactions,LIMO.design.nb_continuous,LIMO.design.zscore,squeeze(LIMO.design.weights(electrode,:,index))',boot_table{electrode});
                    else
                        model = limo_glm1_boot(Y(:,index)',X(index,:),LIMO.design.nb_conditions,LIMO.design.nb_interactions,LIMO.design.nb_continuous,LIMO.design.zscore,squeeze(LIMO.design.weights(electrode,index))',boot_table{electrode});
                    end
                else
                    % index = [1:size(Yr,3)];
                    %                             LIMO.design.weights = squeeze(W(electrode,:));
                    model = limo_glm1_boot(squeeze(Yr(electrode,:,:))',LIMO,boot_table);
                    
                end
                modelTmp{e} = model;
                
                %                                                 H0_Betas(electrode,:,:,:) = model.Betas;
                
            end
            % update the files to be stored on the disk
            for e = 1:size(array,1)
                electrode = array(e);
                model = modelTmp{e};
                
                for B = 1:nboot % now loop because we use cells
                    H0_R2(electrode,:,1,B) = model.R2{B};
                    H0_R2(electrode,:,2,B) = model.F{B};
                    H0_R2(electrode,:,3,B) = model.p{B};
                    
                    if prod(LIMO.design.nb_conditions) ~=0
                        if length(LIMO.design.nb_conditions) == 1
                            tmp_H0_Conditions(electrode,:,1,1,B) = model.conditions.F{B};
                            tmp_H0_Conditions(electrode,:,1,2,B) = model.conditions.p{B};
                        else
                            for i=1:length(LIMO.design.nb_conditions)
                                tmp_H0_Conditions(electrode,:,i,1,B) = model.conditions.F{B}(i,:);
                                tmp_H0_Conditions(electrode,:,i,2,B) = model.conditions.p{B}(i,:);
                            end
                        end
                    end
                    
                    if LIMO.design.fullfactorial == 1
                        if length(LIMO.design.nb_interactions) == 1
                            tmp_H0_Interaction_effect(electrode,:,1,1,:) = model.interactions.F{B};
                            tmp_H0_Interaction_effect(electrode,:,1,2,:) = model.interactions.p{B};
                        else
                            for i=1:length(LIMO.design.nb_interactions)
                                tmp_H0_Interaction_effect(electrode,:,i,1,:) = model.interactions.F{B}(i,:);
                                tmp_H0_Interaction_effect(electrode,:,i,2,:) = model.interactions.p{B}(i,:);
                            end
                        end
                    end
                    
                    if LIMO.design.nb_continuous ~=0
                        if LIMO.design.nb_continuous == 1
                            tmp_H0_Covariates(electrode,:,1,1,B) = model.continuous.F{B};
                            tmp_H0_Covariates(electrode,:,1,2,B) = model.continuous.p{B};
                        else
                            for i=1:LIMO.design.nb_continuous
                                tmp_H0_Covariates(electrode,:,i,1,B) = model.continuous.F{B}(i,:);
                                tmp_H0_Covariates(electrode,:,i,2,B) = model.continuous.p{B}(i,:);
                            end
                        end
                    end
                end
%                 if e == 1 % XXX BEHINGER
%                     cd H0
%                     save boot_table boot_table
%                     save H0_Betas H0_Betas -v7.3
%                     save H0_R2 H0_R2 -v7.3
%                 end
            end
%             close(h)
            warning on;
            
            % save data on the disk and clear out
            cd H0
            save H0_Betas H0_Betas -v7.3
            save H0_R2 H0_R2 -v7.3
            
            
            if prod(LIMO.design.nb_conditions) ~=0
                for i=1:length(LIMO.design.nb_conditions)
                    name = sprintf('H0_Condition_effect_%g',i);
                    H0_Condition_effect = squeeze(tmp_H0_Conditions(:,:,i,:,:));
                    save(name,'H0_Condition_effect','-v7.3');
                    clear H0_Condition_effect
                end
                clear tmp_H0_Conditions
            end
            
            if LIMO.design.fullfactorial == 1
                for i=1:length(LIMO.design.nb_interactions)
                    name = sprintf('H0_Interaction_effect_%g',i);
                    H0_Interaction_effect = squeeze(tmp_H0_Interaction_effect(:,:,i,:,:));
                    save(name,'H0_Interaction_effect','-v7.3');
                    clear H0_Interaction_effect
                end
                clear tmp_H0_Interaction_effect
            end
            
            if LIMO.design.nb_continuous ~=0
                for i=1:length(LIMO.design.nb_continuous)
                    name = sprintf('H0_Covariate_effect_%g',i);
                    H0_Covariate_effect = squeeze(tmp_H0_Covariates(:,:,i,:,:));
                    save(name,'H0_Covariate_effect','-v7.3');
                    clear H0_Covariate_effect
                end
                clear tmp_H0_Covariates
            end
            
            clear electrode model H0_R2; cd ..
            disp(' ');
            
        catch boot_error
            disp('an error occured while attempting to bootstrap the data')
            fprintf('%s \n',boot_error.message);
            rethrow(boot_error)
            return
        end
    end