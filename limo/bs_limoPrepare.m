resp = [];
for k = 1:length(EEG.epoch)
    resp.bsLoc(k) = double(EEG.epoch(k).eventbsLoc{1});
    resp.cNc(k)   = EEG.epoch(k).eventcNc{1};
    resp.trial(k) = EEG.epoch(k).eventtrial{1};
    resp.iStim(k) = EEG.epoch(k).eventiStim{1};
    resp.oStim(k) = EEG.epoch(k).eventoStim{1};
    
end
resp.same =resp.iStim == resp.oStim;

resp.eye(ismember(resp.bsLoc,[1 3])) = 1;%right
resp.eye(ismember(resp.bsLoc,[2 4])) = 2;%left
resp.bs(ismember(resp.bsLoc,[1 4])) = 1; %BS on
resp.bs(ismember(resp.bsLoc,[2 3])) = 0; %BS off

%% selfmade binary designmat
% %% it seems they want the categorical variable in 1 variable...
mat = nan(8,length(resp.bsLoc));

mat(1,:) = resp.eye == 1;
mat(2,:) = resp.eye == 2;
mat(3,:) = resp.bs == 1;
mat(4,:) = resp.bs == 0;
mat(5,:) = resp.cNc == 1;
mat(6,:) = resp.cNc == 2;
resp.spatCont(resp.cNc==2) = ~resp.same(resp.cNc==2);
resp.spatCont(resp.cNc==1) = resp.same(resp.cNc==1);
mat(7,:) = resp.spatCont==1;%resp.same == 1;
mat(8,:) = resp.spatCont==0;%resp.same == 0;
% imagesc(sortrows(mat'))
%% THIS IS WHAT WORKS WITH LIMO 1.4
mat = nan(4,length(resp.bsLoc));
mat(1,:) = resp.eye-1;
mat(2,:) = resp.bs ;
mat(3,:) = ~(resp.cNc-1);
mat(4,:) = resp.spatCont;
testDesignMat = mat;
save('testDesignMat','testDesignMat')

%% plot the selfmade binary designmat
figure
h = imagesc(1-sortrows(mat',[1:8])),colormap('gray')
[row col] = find(diff(sortrows(mat',[1:8]),1,1)~=0)

set(gca,'XtickLabels',{'EyeL','EyeR','iBS','oBS','tempCon','nTempCon','spatCon','nSpatCon'})
k = vline([2.5 4.5 6.5],'r')
l = hline(row,'b')
set([k],'LineWidth',2)
set([k],'LineWidth',4)
% be_print('file','/work/behinger/Dropbox/Masterarbeit/designmat','eps','')
%% plot the non binary designmat
figure
h = imagesc(1-sortrows(mat',[1:4])),colormap('gray')
[row col] = find(diff(sortrows(mat',[1:4]),1,1)~=0)

set(gca,'XtickLabels',{'Eye','BS','tempCon','SpatCon'})
k = vline([2.5 4.5 6.5],'r')
l = hline(row,'b')
set([k],'LineWidth',2)
set([k],'LineWidth',4)
%%

addpath('/home/student/b/behinger/Documents/MATLAB/eeglab_dev/plugins/limo_eeg_v1.3/limo_cluster_functions')
expected_chanlocs = EEG.chanlocs;
channeighbstructmat = tmpNeigh; %limo_neighbourdist(EEG);
save('expected_chanlocs','expected_chanlocs','channeighbstructmat')
