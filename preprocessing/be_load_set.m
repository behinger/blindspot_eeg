function [EEG,p] = be_load_set(p,varargin)
%%  [EEG,p] = be_load_set(p,varargin)
if ischar(p)
    p = be_generate_paths(p);
end
if nargin >= 2
    if ~isnumeric(varargin{1})
        error('Setselection has to be a Number')
    end
    setIdx = varargin{1};
else
    for i = 1:length(p.full.sets)
        fprintf('%i : %s | %s \n',i,p.full.setsDate{i}, p.full.sets{i});
    end
    setIdx=input('Please choose a Set to load: ');
end
if nargin ==3
    loadmode = varargin{2};
else
        loadmode ='all';
end

EEG = pop_loadset('filename',p.full.sets{setIdx},'loadmode',loadmode);


end