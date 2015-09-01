function model = bs_limo_glm1_boot_perm(varargin)

% Boostrapped version of limo_glm1
% Importantly it also runs per electrodes - but do N bootstraps to obtain
% the distributon of F (and associated p values) under H0
% H0 is obtained by either by resampling from centered data (categorical designs)
% or sampling Y but leaving X intact, i.e. breaking the link between Y and X
%
% FORMAT:
% model = limo_glm1_boot(Y,LIMO,boot_table)
% model = limo_glm1_boot(Y,X,nb_conditions,nb_interactions,nb_continuous,zscore,method,boot_table)
%
% INPUTS/OUPUTS: see limo_glm1
%                boot_table is an optional argument - this is the
%                resampling table - if one calls limo_glm1_boot to loop
%                troughout electrodes, this might a good idea to provide
%                such table so that the same resampling applies to each
%                electrodes
%
% See also
% LIMO_DESIGN_MATRIX, LIMO_WLS, LIMO_IRLS, LIMO_EEG(4)
%
% Cyril Pernet v1 18-07-2012
% -----------------------------
%  Copyright (C) LIMO Team 2012

%% varagin
nboot = 599; %

if nargin == 2 || nargin == 3
    y               = varargin{1};
    X               = varargin{2}.design.X;
    nb_conditions   = varargin{2}.design.nb_conditions;
    nb_interactions = varargin{2}.design.nb_interactions;
    nb_continuous   = varargin{2}.design.nb_continuous;
    z               = varargin{2}.design.zscore;
    method          = varargin{2}.design.method;
    if nargin == 2
        boot_table = randi(size(Y,1),size(Y,1),nboot);
    elseif nargin == 3
        boot_table = varargin{3};
        nboot = size(boot_table,2);
    end
elseif nargin == 7 || nargin == 8
    y               = varargin{1};
    X               = varargin{2};
    nb_conditions   = varargin{3};
    nb_interactions = varargin{4};
    nb_continuous   = varargin{5};
    z               = varargin{6};
    method          = varargin{7};
    if nargin == 7
        boot_table = randi(size(Y,1),size(Y,1),nboot);
    elseif nargin == 8
        boot_table = varargin{8};
        nboot = size(boot_table,2);
    end
else
    error('varargin error in limo_glm1_boot')
end

clear varargin
nb_factors = numel(nb_conditions);
if nb_factors == 1 && nb_conditions == 0
    nb_factors = 0;
end

% -----------
%% Data check
% -----------

if size(y,1)~=size(X,1)
    error('The number of events in Y and the design matrix are different')
end

if nb_interactions == 0
    nb_interactions = [];
end

% ----------
%% Bootstrap
% -----------
design = X;

% if categorical design, center data 1st
% ---------------------------------------
% if nb_continuous == 0
%     centered_y = NaN(size(y,1),size(y,2));
%     if ~isempty(nb_interactions) && length(nb_interactions)<length(nb_conditions)
%         [tmpX interactions] = limo_make_interactions(X(:,1:(end-1-length(nb_interactions))), nb_conditions);
%         if length(interactions) == 1
%             start_at = sum(nb_conditions);
%         else
%             start_at = sum(nb_conditions)+sum(interactions(1:end-1));
%         end
%         
%         for cel=(start_at+1):(start_at+interactions(end))
%             index = find(tmpX(:,cel));
%             centered_y(index,:) = y(index,:) - repmat(mean(y(index,:),1),[size(y(index,:),1)],1);
%         end
%         
%         
%     elseif ~isempty(nb_interactions)
%         
%         % look up the last interaction to get unique groups
%         if length(nb_interactions) == 1
%             start_at = sum(nb_conditions);
%         else
%             start_at = sum(nb_conditions)+sum(nb_interactions(1:end-1));
%         end
%         
%         for cel=(start_at+1):(start_at+nb_interactions(end))
%             index = find(X(:,cel));
%             centered_y(index,:) = y(index,:) - repmat(mean(y(index,:),1),length(index),1);
%         end
%         
%     elseif size(nb_conditions,2) == 1
%         % no interactions because just 1 factor
%         for cel=1:nb_conditions
%             index = find(X(:,cel));
%             centered_y(index,:) = y(index,:) - repmat(mean(y(index,:),1),length(index),1);
%         end
%         
%     else
%         % create fake interaction to get groups
%         [tmpX interactions] = limo_make_interactions(X(:,1:(end-1)), nb_conditions);
%         if length(interactions) == 1
%             start_at = sum(nb_conditions);
%         else
%             start_at = sum(nb_conditions)+sum(interactions(1:end-1));
%         end
%         
%         for cel=(start_at+1):(start_at+interactions(end))
%             index = find(tmpX(:,cel));
%             centered_y(index,:) = y(index,:) - repmat(mean(y(index,:),1),[size(y(index,:),1)],1);
%         end
%     end
%     clear y
% end

% compute for each bootstrap
% ---------------------------
for B = 1:nboot
    
    % create data under H0
%     if nb_continuous == 0
%         % sample from the centered data in categorical designs
%         Y = centered_y(boot_table(:,B),:);
%         X = design(boot_table(:,B),:); % resample X as Y
        
%     else
        % sample and break the link between Y and (regression and AnCOVA designs)
        Y = y(boot_table(:,B),:);
        X = design;
        if z == 1 % rezscore the covariates
            N = nb_conditions + nb_interactions;
            if N==0
                if sum(mean(X(:,1:end-1),1)) > 10e-15
                    X(:,1:end-1) = zscore(X(:,1:end-1));
                end
            else
                if sum(mean(X(:,N+1:end-1),1)) > 10e-15
                    X(:,N+1:end-1) = zscore(X(:,N+1:end-1));
                end
            end
        end
%     end
%     
    % ------------------------------
    % Compute F for dummy variables
    % ------------------------------
    
    % -------------------------
    if nb_factors == 1   %  1-way ANOVA
        % -------------------------
        
        % total sum of squares, projection matrix for errors, residuals and betas
        % -----------------------------------------------------------------------
        T     = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));  % SS Total
        R     = eye(size(Y,1)) - (X*pinv(X));                                      % Projection on E
        E     = (Y'*R*Y);                                                          % SS Error
        
        % compute Beta parameters
        if strcmp(method,'OLS')
            Betas = pinv(X)*Y;
        elseif strcmp(method,'WLS')
            [Betas,W] = limo_WLS(X,Y);
        elseif strcmp(method,'IRLS')
            [Betas,W] = limo_IRLS(X,Y);
        end
        model.Betas(:,:,B) = Betas';
        
        % compute model R^2
        % -----------------
        C = eye(size(X,2));
        C(:,size(X,2)) = 0;
        C0 = eye(size(X,2)) - C*pinv(C);
        X0 = X*C0;  % Reduced model
        R0 = eye(size(Y,1)) - (X0*pinv(X0));
        M  = R0 - R;  % Projection matrix onto Xc
        H  = (Betas'*X'*M*X*Betas);  % SS Effect
        Rsquare   = diag(H)./diag(T); % Variances explained
        F_Rsquare = (diag(H)./(rank(X)-1)) ./ (diag(E)/(size(Y,1)-rank(X)));
        p_Rsquare = 1 - fcdf(F_Rsquare, (rank(X)-1), (size(Y,1)-rank(X)));
        
        % compute F for categorical variables
        % -----------------------------------
        if nb_conditions ~= 0 && nb_continuous == 0
            df_conditions   = rank(C)-1;
            F_conditions    = F_Rsquare;
            pval_conditions = p_Rsquare;
            
        elseif nb_conditions ~= 0 && nb_continuous ~= 0
            C = eye(size(X,2));
            C(:,(nb_conditions+1):size(X,2)) = 0;
            C0 = eye(size(X,2)) - C*pinv(C);
            X0 = X*C0; % Here the reduced model includes the covariates
            R0 = eye(size(Y,1)) - (X0*pinv(X0));
            M  = R0 - R;
            H  = (Betas'*X'*M*X*Betas);
            df_conditions = rank(C)-1;
            F_conditions    = (diag(H)/(rank(C)-1)) ./ (diag(E)/(size(Y,1)-rank(X)));
            pval_conditions = 1 - fcdf(F_conditions(:), df_conditions, (size(Y,1)-rank(X)));
        end
        
        model.conditions.F{B}  = F_conditions;
        model.conditions.p{B}  = pval_conditions;
        
        % ------------------------------------------------
    elseif nb_factors > 1  && isempty(nb_interactions) % N-ways ANOVA without interactions
        % ------------------------------------------------
        
        % compute basic SS total, projection matrices and parameters
        T        = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));
        R        = eye(size(Y,1)) - (X*pinv(X));
        E        = (Y'*R*Y);
        % compute Beta parameters with weights
        if strcmp(method,'OLS')
            Betas = pinv(X)*Y;
        elseif strcmp(method,'WLS')
            [Betas,W] = limo_WLS(X,Y);
        elseif strcmp(method,'IRLS')
            [Betas,W] = limo_IRLS(X,Y);
        end
        model.Betas(:,:,B) = Betas';
        
        % --------------------
        % compute model R^2
        % --------------------
        C = eye(size(X,2));
        C(:,size(X,2)) = 0;
        C0   = eye(size(X,2)) - C*pinv(C);
        X0   = X*C0; % Reduced model (i.e. only intercept)
        R0   = eye(size(Y,1)) - (X0*pinv(X0));
        M    = R0 - R;      % M is the projection matrix onto Xc
        H    = (Betas'*X'*M*X*Betas);   % SSCP Hypothesis (Effect)
        Rsquare   = diag(H)./diag(T); % Variances explained per Y
        F_Rsquare = (diag(H)./(rank(X)-1)) ./ (diag(E)/(size(Y,1)-rank(X)));
        p_Rsquare = 1 - fcdf(F_Rsquare, (rank(X)-1), (size(Y,1)-rank(X)));
        
        % --------------------------------------
        % compute F and p values of each factor
        % --------------------------------------
        
        df_conditions = zeros(1,length(nb_conditions));
        F_conditions = zeros(length(nb_conditions),size(Y,2));
        pval_conditions = zeros(length(nb_conditions),size(Y,2));
        
        eoi = zeros(1,size(X,2));
        eoi(1:nb_conditions(1)) = 1:nb_conditions(1);
        eoni = [1:size(X,2)];
        eoni = find(eoni - eoi);
        
        for f = 1:length(nb_conditions)
            C = eye(size(X,2));
            C(:,eoni) = 0;
            C0   = eye(size(X,2)) - C*pinv(C);
            X0   = X*C0;
            R0   = eye(size(Y,1)) - (X0*pinv(X0));
            M    = R0 - R;
            H    = (Betas'*X'*M*X*Betas);
            df_conditions(f) = rank(C)-1;
            F_conditions(f,:)    = (diag(H)/df_conditions(f)) ./ (diag(E)/(size(Y,1)-rank(X)));
            pval_conditions(f,:) = 1 - fcdf(F_conditions(f,:), df_conditions(f), (size(Y,1)-rank(X)));
            
            % update factors
            if f<length(nb_conditions)
                update = max(find(eoi));
                eoi = zeros(1,size(X,2));
                eoi((update+1):(update+nb_conditions(f+1))) = update + (1:nb_conditions(f+1));
                eoni = [1:size(X,2)];
                eoni = find(eoni - eoi);
            end
        end
        model.conditions.F{B}  = F_conditions;
        model.conditions.p{B}  = pval_conditions;
        
        
        % ------------------------------------------------
    elseif nb_factors > 1  && ~isempty(nb_interactions) % N-ways ANOVA with interactions
        % ------------------------------------------------
        
        % compute basic SS total, projection matrices and parameters
        T        = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));
        R        = eye(size(Y,1)) - (X*pinv(X));
        E        = (Y'*R*Y);
        % compute Beta parameters with weights
        if strcmp(method,'OLS')
            Betas = pinv(X)*Y;
        elseif strcmp(method,'WLS')
            [Betas,W] = limo_WLS(X,Y);
        elseif strcmp(method,'IRLS')
            [Betas,W] = limo_IRLS(X,Y);
        end
        model.Betas(:,:,B) = Betas';
        
        % --------------------
        % compute model R^2
        % --------------------
        C = eye(size(X,2));
        C(:,size(X,2)) = 0;
        C0   = eye(size(X,2)) - C*pinv(C);
        X0   = X*C0; % Reduced model (i.e. only intercept)
        R0   = eye(size(Y,1)) - (X0*pinv(X0));
        M    = R0 - R;      % M is the projection matrix onto Xc
        H    = (Betas'*X'*M*X*Betas);   % SSCP Hypothesis (Effect)
        Rsquare   = diag(H)./diag(T); % Variances explained per Y
        F_Rsquare = (diag(H)./(rank(X)-1)) ./ (diag(E)/(size(Y,1)-rank(X)));
        p_Rsquare = 1 - fcdf(F_Rsquare, (rank(X)-1), (size(Y,1)-rank(X)));
        
        
        % ---------------------------------------------------
        % start by ANOVA without interaction for main effects
        % ---------------------------------------------------
        
        df_conditions = zeros(1,length(nb_conditions));
        F_conditions = zeros(length(nb_conditions),size(Y,2));
        pval_conditions = zeros(length(nb_conditions),size(Y,2));
        
        % covariates
        covariate_columns = [(sum(nb_conditions)+sum(nb_interactions)+1):(size(X,2)-1)];
        
        % main effects
        dummy_columns = 1:sum(nb_conditions);
        
        % re-define X
        x = [X(:,dummy_columns) X(:,covariate_columns) ones(size(X,1),1)];
        
        % run same model as above
        R        = eye(size(Y,1)) - (x*pinv(x));
        % compute Beta parameters with weights
        if strcmp(method,'OLS')
            betas = pinv(x)*Y;
        elseif strcmp(method,'WLS')
            [betas,W] = limo_WLS(x,Y);
        elseif strcmp(method,'IRLS')
            [betas,W] = limo_IRLS(x,Y);
        end
        
        eoi = zeros(1,size(x,2));
        eoi(1:nb_conditions(1)) = 1:nb_conditions(1);
        eoni = [1:size(x,2)];
        eoni = find(eoni - eoi);
        
        for f = 1:length(nb_conditions)
            C = eye(size(x,2));
            C(:,eoni) = 0;
            C0   = eye(size(x,2)) - C*pinv(C);
            X0   = x*C0;
            R0   = eye(size(Y,1)) - (X0*pinv(X0));
            M    = R0 - R;
            H(f,:) = diag((betas'*x'*M*x*betas));
            df_conditions(f) = rank(C)-1;
            F_conditions(f,:)    = (H(f,:)./df_conditions(f)) ./ (diag(E)./(size(Y,1)-rank(X)))';
            pval_conditions(f,:) = 1 - fcdf(F_conditions(f,:), df_conditions(f), (size(Y,1)-rank(X)));
            
            % update factors
            if f<length(nb_conditions)
                update = max(find(eoi));
                eoi = zeros(1,size(x,2));
                eoi((update+1):(update+nb_conditions(f+1))) = update + (1:nb_conditions(f+1));
                eoni = [1:size(x,2)];
                eoni = find(eoni - eoi);
            end
        end
        model.conditions.F{B}  = F_conditions;
        model.conditions.p{B}  = pval_conditions;
        
        % ---------------------------
        % now deal with interactions
        % ---------------------------
        
        if nb_factors == 2 && nb_continuous == 0 % the quick way with only one interaction
            HI = diag(T)' - H(1,:) - H(2,:) - diag(E)';
            df_interactions = prod(df_conditions);
            F_interactions  = (HI./df_interactions) ./ (diag(E)/(size(Y,1)-rank(X)))';
            pval_interactions  = 1 - fcdf(F_interactions, df_interactions, (size(Y,1)-rank(X)));
            
        else % run through each interaction
            
            % part of X unchanged
            Main_effects = [X(:,dummy_columns)];
            Cov_and_Mean = [X(:,covariate_columns) ones(size(Y,1),1)];
            
            % get interactions
            start = size(Main_effects,2)+1;
            for i=1:length(nb_interactions)
                I{i} = X(:,start:(start+nb_interactions(i)-1));
                start = start+nb_interactions(i);
            end
            start = size(Main_effects,2)+1;
            
            % check interaction levels
            index = 1;
            for n=2:nb_factors
                combinations = nchoosek([1:nb_factors],n); % note it matches I above because computed with nchoosek the same way in limo_design_matrix
                for c = 1:size(combinations,1)
                    interaction{index} = combinations(c,:);
                    index = index + 1;
                end
            end
            
            add = 0; start_at_I = 1;
            % run substituting and/or incrementing parts of X
            for f = 1:length(nb_interactions)
                
                % re-define X with interactions
                test = size(interaction{f},2);
                if test == 2
                    x = [Main_effects I{f} Cov_and_Mean];
                    add = add+1;
                else
                    if add == test
                        for a = start_at_I:add
                            Main_effects = [Main_effects I{a}];
                        end
                        start = size(Main_effects,2)+1;
                        start_at_I = add+1;
                    end
                    x = [Main_effects I{f} Cov_and_Mean];
                end
                
                % run same model as above
                R  = eye(size(Y,1)) - (x*pinv(x));
                if strcmp(method,'OLS')
                    betas = pinv(x)*Y;
                elseif strcmp(method,'WLS')
                    [betas,W] = limo_WLS(x,Y);
                elseif strcmp(method,'IRLS')
                    [betas,W] = limo_IRLS(x,Y);
                end
                
                
                eoi = zeros(1,size(x,2));
                eoi(start:(start-1+nb_interactions(f))) = start:(start-1+nb_interactions(f));
                eoni = [1:size(x,2)];
                eoni = find(eoni - eoi);
                
                C = eye(size(x,2));
                C(:,eoni) = 0;
                C0   = eye(size(x,2)) - C*pinv(C);
                X0   = x*C0;
                R0   = eye(size(Y,1)) - (X0*pinv(X0));
                M    = R0 - R;
                HI(f,:) = diag((betas'*x'*M*x*betas))';
            end
            
            % get appropriate df and F/p values
            df_interactions = zeros(1,length(nb_interactions));
            F_interactions = zeros(length(nb_interactions),size(Y,2));
            pval_interactions = zeros(length(nb_interactions),size(Y,2));
            
            for f = 1:length(nb_interactions)
                dfs = df_conditions(interaction{f});
                df_interactions(f) = prod(dfs);
                F_interactions(f,:) = (HI(f,:)./df_interactions(f)) ./ (diag(E)/(size(Y,1)-rank(X)))';
                pval_interactions(f,:) = 1 - fcdf(F_interactions(f,:), df_interactions(f), (size(Y,1)-rank(X)));
            end
        end
        model.interactions.F{B}  = F_interactions;
        model.interactions.p{B}  = pval_interactions;
    end
    
    
    % -----------------------------------
    %% compute F for continuous variables
    % -----------------------------------
    
    if nb_continuous ~=0
        
        if nb_factors == 0
            T     = (Y-repmat(mean(Y),size(Y,1),1))'*(Y-repmat(mean(Y),size(Y,1),1));
            R     = eye(size(Y,1)) - (X*pinv(X));
            E     = (Y'*R*Y);
            % compute Beta parameters with weights Y was resampled, so
            % resample W as well - however X unchanged so use W
            if strcmp(method,'OLS')
                Betas = pinv(X)*Y;
            elseif strcmp(method,'WLS')
                [Betas,W] = limo_WLS(X,Y);
            elseif strcmp(method,'IRLS')
                [Betas,W] = limo_IRLS(X,Y);
            end
            model.Betas(:,:,B) = Betas';
            
            % compute model R^2
            % -----------------
            C = eye(size(X,2));
            C(:,size(X,2)) = 0;
            C0 = eye(size(X,2)) - C*pinv(C);
            X0 = X*C0;
            R0 = eye(size(Y,1)) - (X0*pinv(X0));
            M  = R0 - R;
            H  = (Betas'*X'*M*X*Betas);
            Rsquare   = diag(H)./diag(T);
            F_Rsquare = (diag(H)./(rank(X)-1)) ./ (diag(E)/(size(Y,1)-rank(X)));
            p_Rsquare = 1 - fcdf(F_Rsquare, (rank(X)-1), (size(Y,1)-rank(X)));
        end
        
        if nb_factors == 0 && nb_continuous == 1 % simple regression
            model.continuous.F{B}  = F_Rsquare;
            model.continuous.p{B}  = p_Rsquare;
            
        else  % ANCOVA
            
            % pre-allocate space
            F_continuous = zeros(nb_continuous,size(Y,2));
            pval_continuous = zeros(nb_continuous,size(Y,2));
            
            % compute
            N_conditions = sum(nb_conditions) + sum(nb_interactions);
            for n = 1:nb_continuous
                C    = zeros(size(X,2));
                C(N_conditions+n,N_conditions+n) = 1;
                C0   = eye(size(X,2)) - C*pinv(C);
                X0   = X*C0;
                R0   = eye(size(Y,1)) - (X0*pinv(X0));
                M    = R0 - R;
                H    = Betas'*X'*M*X*Betas;
                F_continuous(n,:) = (diag(H)./(rank(C))) ./ (diag(E)/(size(Y,1)-rank(X)));
                pval_continuous(n,:) = 1 - fcdf(F_continuous(n,:), 1, (size(Y,1)-rank(X)));
            end
            model.continuous.F{B}  = F_continuous';
            model.continuous.p{B}  = pval_continuous';
        end
    end
    
    % ----------------------------
    %% update the model structure
    % ----------------------------
    model.R2{B} = Rsquare;
    model.F{B} = F_Rsquare;
    model.p{B} = p_Rsquare;
end
