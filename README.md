# blindspot_eeg
Published in JoN 2015

Disclaimer: These scripts are unorganized, incomplete and not as good documented as they could/should be

I put them online anyway because 
  a) I'm interested whether I made any mistakes that can influence my conclusions I drew from the data
  b) It might be interesting for other scientists/data analysts to see how we used TFCE/LIMO
  
I am very open for suggestions, comments and questions.


Here is a short description of the analysis pipeline, followed by a description of the LIMO statistical stages.


==Analysis==
>> EEG = be_import(p); % Imports the file using pop_import

>> EEG = bs_add_events_all_v2(EEG,p); % Adds the events from the eye tracker

>> EEG = be_resample(EEG,p,500); % downsample to 500 Hz

>> EEG = be_load_xensor(EEG,p); % Add the xensor individual electrode position

>> EEG = be_reject_channel(EEG,p,1); % Reject previously (manually) marked channels

>> EEG = be_reref(EEG,p); % Average Ref AFTER channel rejection

>> EEG = be_load_ICA(EEG,p); % At this stage I have the ICA already run,
I usually do a 1 Hz Highpass for the ICA, let it run, extract the
weights and put them on the unfiltered data (do the 1Hz even so you will
not need to filter in the first place)

>> EEG = be_filter_cont(EEG,p,0.1,120); % bandpass 0.1 to 120, 0.1 could
also be left out or put up to 1Hz depending on analysis.

>> EEG = be_clean_continuous(EEG,p,1); % Here my previous automatic
cleaning times are loaded and the bad data is rejected. I clean the data
manually

>> EEG = be_ICA_mark(EEG,p,1); % My ICA bad components a) identified by
eyetracking, b) by automatic precedures. For eyetracking combined data some Algorithm works very well (Plöchl, Ossandon, König 201?, see also Olaf Dimigens Eyetracking/EEG toolbox).

>> EEG = pop_subcomp(EEG,[],0); % I delete the previously marked components

>> EEG = be_epoch(EEG,p,'etFixOnset',[-1.5 1],[-300 -100]); % The data
are epoched and baseline corrected

>> EEG = pop_interp(EEG,chanAll,'spherical'); % Interpolate missing channels


==LIMO==
t.b.d.
