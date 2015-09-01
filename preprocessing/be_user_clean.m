k = 20;
p = flags(k );
EEG = be_import(p);
EEG = be_resample(EEG,p,500);
EEG = be_filter_cont(EEG,p,1,120);
% EEG = be_reject_channel(EEG,p);
% EEG = be_reref(EEG,p);
be_clean_continuous(EEG,p);
