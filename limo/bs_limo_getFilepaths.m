function [folderPath,fileNames,fullPath] = bs_limo_getFilepaths(varargin)

g = be_inputcheck(varargin, ...
{'filename','string','',''; 
    'subjects','integer',[1:10],1:10;
    'path','string','','/net/store/nbp/EEG/blind_spot/data/limoTest'});
if ischar(g)
    error(g)
end
if ~exist(g.path, 'dir')
    error('%s does not exist',g.path)
end

subjList = {'A','B','C','D','E','G','H','J','K','L','M','N','T','U','V'};
for k = 1:length(g.subjects)
   folderPath{k} = [g.path filesep subjList{k}];
   fileNames{k} = g.filename;
   fullPath{k} = [folderPath{k} filesep fileNames{k}];
end

% tmpStr =  sprintf('''/net/store/nbp/EEG/blind_spot/data/limoTest/%s'',',) ;
%  eval(['{' tmpStr(1:end-1) '}']);