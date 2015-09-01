function LIMO = bs_limo_designMat(LIMO,varargin)
%% limo = bs_limo_designMat(limo)
% Input: the limo structure, EEGfile in LIMO.data.data/data_dir is loaded
% and the designmatrix is created.


g = be_inputcheck(varargin, ...
    {'bs','string',{'yes','no'},'yes';
    'session','string',{'yes','no'},'no';
    'interaction','string',{'man','yes','no','full'},'no'; %overwrites everthing
    'full','string',{'yes','no'},'no';
    'stim','integer',[1 2],2});
if ischar(g)
    error(g)
end
EEG = pop_loadset('filename', LIMO.data.data,'filepath',LIMO.data.data_dir);

resp = [];

for k = 1:length(EEG.epoch)
    if isempty(EEG.epoch(k).eventbsLoc{1})
        resp.bsLoc(k) = nan;
        resp.cNc(k) = nan;
        resp.trial(k) = nan;
        resp.iStim(k) = nan;
        resp.oStim(k) = nan;
        
        continue
    end
    resp.bsLoc(k) = double(EEG.epoch(k).eventbsLoc{1});
    resp.cNc(k)   = EEG.epoch(k).eventcNc{1};
    resp.trial(k) = EEG.epoch(k).eventtrial{1};
    resp.iStim(k) = EEG.epoch(k).eventiStim{1};
    resp.oStim(k) = EEG.epoch(k).eventoStim{1};
    
end
resp.spatCont = nan(1,length(resp.bsLoc));
resp.leftRight = nan(1,length(resp.bsLoc));
resp.bs = nan(1,length(resp.bsLoc));
resp.same =resp.iStim == resp.oStim;
if g.stim == 1
    resp.spatCont = ~resp.same;
else
    resp.spatCont(resp.cNc==2) = resp.same(resp.cNc==2);
    resp.spatCont(resp.cNc==1) = ~resp.same(resp.cNc==1);
end
resp.leftRight(ismember(resp.bsLoc,[1 3])) = 1;%left
resp.leftRight(ismember(resp.bsLoc,[2 4])) = 2;%right
resp.bs(ismember(resp.bsLoc,[1 4])) = 1; %BS on
resp.bs(ismember(resp.bsLoc,[2 3])) = 0; %BS off

if strcmp(g.session,'yes')
    idxTmp = find(~isnan(resp.trial));
    
    sesIdx = find(abs(diff(resp.trial(idxTmp)))>1000)+1;
    sesIdx = idxTmp(sesIdx);
end

notrial = []; %which to mark with nans?
mat(1,:) = resp.leftRight-1;     % left == 0 right == 1
LIMO.design.XDesc{1} = 'Left=0,right=1';



if strcmp(g.full,'yes')
    mat(end+1,:) = resp.bs ;       % bs, 0 no bs
    LIMO.design.XDesc{end+1} = 'oBs=0,wBs=1';
    if g.stim == 2 %|| ~isempty(strfind(LIMO.data.data_dir,'saccade'))
        mat(end+1,:) = resp.cNc-1;  %temporal contrast, 0 no contrast
        LIMO.design.XDesc{end+1} = 'NoTempCont=0,wTempCont=1';
    end
    mat(end+1,:) = resp.spatCont;  %spatial contrast, 0 no contrast
    LIMO.design.XDesc{end+1} = 'NoSpatCont=0,wSpatCont=1';
    
    
    
else
    if strcmp(g.bs,'yes')
        mat(end+1,:) = resp.bs ;       % bs, 0 no bs
        LIMO.design.XDesc{end+1} = 'oBs=0,wBs=1';
        if g.stim ~=1
            notrial = [notrial find(resp.cNc-1 ~= resp.spatCont)];
        end
    else
        mat(end+1,:) = resp.spatCont;  %spatial contrast, 0 no contrast
        LIMO.design.XDesc{end+1} = 'NoSpatCont=0,wSpatCont=1';
        notrial = [notrial find(resp.bs)];
    end
    if g.stim ==1
        %     mat(end+1,:) = resp.
        %     mat(end+1,:) = resp.spatCont;  %spatial contrast, 0 no contrast
        %     LIMO.design.XDesc{end+1} = 'NoSpatCont=0,wSpatCont=1';
%         notrial = [notrial find(resp.same==0)]; %throw out all, that are not the same
    else
        mat(end+1,:) = resp.cNc-1;  %temporal contrast, 0 no contrast
        LIMO.design.XDesc{end+1} = 'NoTempCont=0,wTempCont=1';
    end
end
if strcmp(g.session,'yes')
    resp.session = [zeros(1,sesIdx) ones(1,length(resp.trial)-sesIdx)];
    mat(end+1,:) = resp.session ;       % bs, 0 no bs
    LIMO.design.XDesc{end+1} = 'session1=0,session2=1';
end

mat(:,notrial) = nan;
if strcmp(g.interaction,'yes') || strcmp(g.interaction,'no') || strcmp(g.interaction,'man')
    LIMO.design.fullfactorial = 0;
elseif strcmp(g.interaction,'full')
    LIMO.design.fullfactorial = 1;
end
LIMO.data.Cat = mat;
[LIMO.design.X, LIMO.design.nb_conditions, LIMO.design.nb_interactions, LIMO.design.nb_continuous] = ...
    limo_design_matrix(EEG.data, LIMO, 0);

LIMO.design.name = 'blindspotGLM';
% interaction
%BS Interaction
if strcmp(g.interaction,'man') || strcmp(g.interaction,'yes')
    
    dispLeft = 2*find(strcmp('Left=0,right=1',LIMO.design.XDesc))-1;
    dispRight = dispLeft+1;
    bsOff = 2*find(strcmp('oBs=0,wBs=1',LIMO.design.XDesc))-1;
    bsOn = bsOff+1;
    tcOff = 2*find(strcmp('NoTempCont=0,wTempCont=1',LIMO.design.XDesc))-1;
    tcOn = tcOff+1;
    scOff = 2*find(strcmp('NoSpatCont=0,wSpatCont=1',LIMO.design.XDesc))-1;
    scOn = scOff+1;
    mat = LIMO.design.X;
    LIMO.design.nb_interactions = [];
    sizeBefore = size(mat,2);
    if ~isempty(bsOn) &&  ~isempty(tcOn) % BS Tempo Interact
        mat = [mat(:,1:end-1) (mat(:,bsOff).*mat(:,tcOff)) mat(:,end)] ;
        LIMO.design.XDesc = [LIMO.design.XDesc {'oBS/noTempospatio Int'}];
        mat = [mat(:,1:end-1) (mat(:,bsOff).*mat(:,tcOn)) mat(:,end)] ;
        LIMO.design.XDesc = [LIMO.design.XDesc {'oBS/Tempospatio Int'}];
        mat = [mat(:,1:end-1) (mat(:,bsOn).*mat(:,tcOff)) mat(:,end)] ;
        LIMO.design.XDesc = [LIMO.design.XDesc {'iBS/noTempospatio'}];
        mat = [mat(:,1:end-1) (mat(:,bsOn).*mat(:,tcOn)) mat(:,end)] ;
        LIMO.design.XDesc = [LIMO.design.XDesc {'iBS/Tempospatio'}];
        LIMO.design.nb_interactions = [LIMO.design.nb_interactions 4];
    end
    if ~isempty(scOff) && g.stim ~=1 % Temp / Spatio Interact
        mat = [mat(:,1:end-1) (mat(:,tcOff).*mat(:,scOff)) mat(:,end)] ;
        LIMO.design.XDesc = [LIMO.design.XDesc {'noTemp/noSpat Int'}];
        mat = [mat(:,1:end-1) (mat(:,tcOff).*mat(:,scOn)) mat(:,end)] ;
        LIMO.design.XDesc = [LIMO.design.XDesc {'noTemp/Spat Int'}];
        mat = [mat(:,1:end-1) (mat(:,tcOn).*mat(:,scOff)) mat(:,end)] ;
        LIMO.design.XDesc = [LIMO.design.XDesc {'Temp/noSpat Int'}];
        mat = [mat(:,1:end-1) (mat(:,tcOn).*mat(:,scOn)) mat(:,end)] ;
        LIMO.design.XDesc = [LIMO.design.XDesc {'Temp/Spat Int'}];
        LIMO.design.nb_interactions = [LIMO.design.nb_interactions 4];
        
    end
    if ~isempty(bsOn)  &&  ~isempty(scOn) && g.stim ~=1% BS / Spatio Interact
        mat = [mat(:,1:end-1) (mat(:,bsOff).*mat(:,scOff)) mat(:,end)] ;
        LIMO.design.XDesc = [LIMO.design.XDesc {'oBS/noSpatio Int'}];
        mat = [mat(:,1:end-1) (mat(:,bsOff).*mat(:,scOn)) mat(:,end)] ;
        LIMO.design.XDesc = [LIMO.design.XDesc {'oBS/Spatio Int'}];
        mat = [mat(:,1:end-1) (mat(:,bsOn).*mat(:,scOff)) mat(:,end)] ;
        LIMO.design.XDesc = [LIMO.design.XDesc {'iBS/noSpatio'}];
        mat = [mat(:,1:end-1) (mat(:,bsOn).*mat(:,scOn)) mat(:,end)] ;
        LIMO.design.XDesc = [LIMO.design.XDesc {'iBS/Spatio'}];
        LIMO.design.nb_interactions = [LIMO.design.nb_interactions 4];
    end
    if g.stim ==1 % Disp/BS Interact
        mat = [mat(:,1:end-1) (mat(:,bsOff).*mat(:,dispLeft)) mat(:,end)] ;
        LIMO.design.XDesc = [LIMO.design.XDesc {'oBS/leftDisp Int'}];
        mat = [mat(:,1:end-1) (mat(:,bsOff).*mat(:,dispRight)) mat(:,end)] ;
        LIMO.design.XDesc = [LIMO.design.XDesc {'oBS/rightDisp Int'}];
        mat = [mat(:,1:end-1) (mat(:,bsOn).*mat(:,dispLeft)) mat(:,end)] ;
        LIMO.design.XDesc = [LIMO.design.XDesc {'iBS/leftDisp'}];
        mat = [mat(:,1:end-1) (mat(:,bsOn).*mat(:,dispRight)) mat(:,end)] ;
        LIMO.design.XDesc = [LIMO.design.XDesc {'iBS/rightDisp'}];
        LIMO.design.nb_interactions = [LIMO.design.nb_interactions 4];
    end
    
    LIMO.design.X = mat;
    %     LIMO.design.nb_conditions = [LIMO.design.nb_conditions ones(1,size(mat,2)-sizeBefore)];
    
    
    %     LIMO.design.nb_conditions = [LIMO.design.nb_conditions];
    Betas = zeros(size(LIMO.data.chanlocs,2),size(EEG.data,2),size(LIMO.design.X,2)); save Betas Betas; clear Betas
end

% LIMO = bs_limo_designMat_removeDouble(LIMO);

