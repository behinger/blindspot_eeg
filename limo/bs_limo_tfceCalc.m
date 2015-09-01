function [mask M]  = bs_limo_tfceCalc(MCC,FileName,recalcMCC,LIMO,toplot)

if MCC ~= 5
    return
end
mask = [];
M = [];
cd([LIMO.dir '/groupLimo'])
if exist(sprintf('H0/mcc_%i_%s',MCC,FileName),'file') && ~recalcMCC
%     load(sprintf('H0/mcc_%i_%s',MCC,FileName))
elseif MCC == 5
    tmpIdx = regexp(FileName,'(_)');
    [~,tfce_name] = fileparts(['tfce_' FileName(1:tmpIdx(2)-1)]);
    [~,tfce_H0_name] = fileparts(['tfce_H0_'  FileName(1:tmpIdx(2)-1)]);
    [~,tfce_nameFull] = fileparts(['tfce_' FileName]);
    [~,tfce_H0_nameFull] = fileparts(['tfce_H0_'  FileName]);
    %         tfce_name = sprintf('tfce_paired_samples_ttest_parameter_%s',num2str(pred')');
    %         tfce_H0_name = sprintf('tfce_H0_paired_samples_ttest_parameter_%s',num2str(pred')');
    if ~exist('tfce','dir')
        mkdir('tfce')
    end
    
    
    fprintf('Thresholding  Sample  using TFCE \n');
    tmp = limo_tfce(squeeze(toplot),LIMO.data.neighbouring_matrix);
    eval([tfce_name '=tmp;']) % sorry, its the LIMO toolbox way...
    save(['tfce' filesep tfce_nameFull],tfce_name, '-v7.3'); clear('tmp',tfce_name);
    % do tfce for the boostrapped data
    fprintf('Thresholding H0 Sample using TFCE \n');
    tmpH0Dat = load(['H0/H0_' FileName]);
    fn = fieldnames(tmpH0Dat);
    if length(fn) ~= 1, error('The Statisticsfile should only contain one matrix: %s',['H0/H0_' FileName]);end;
    tmp = limo_tfce(squeeze(tmpH0Dat.(fn{1})(:,:,1,:)),LIMO.data.neighbouring_matrix);
    eval([tfce_H0_name '=tmp;']) % sorry, its the LIMO toolbox way...
    tmpH0Dat =[];
    save(['H0', filesep, tfce_H0_nameFull],tfce_H0_name, '-v7.3'); clear('tmp',tfce_H0_name);
end
