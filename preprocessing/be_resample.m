function [EEG] = be_resample(EEG,p,varargin)
    if nargin < 3
        newSamp = 512;
    elseif nargin == 3
        newSamp = varargin{1};
    else
        error('Nargin seems to be wrong %i',nargin)
    end
        
p = be_generate_paths(p);
if ~check_EEG(EEG.preprocess,'Resample')
    tmpSetname = EEG.setname;
    EEG = pop_resample( EEG, newSamp);
    EEG.setname=tmpSetname;
    EEG.preprocessInfo.resampleDate = datestr(now);
    EEG.preprocess = [EEG.preprocess 'Resample'];
end