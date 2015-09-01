# blindspot_eeg
Published in JoN 2015

Predictions of Visual Content across Eye Movements and Their Modulation by Inferred Information

    Benedikt V. Ehinger,  Peter König and José P. Ossandón
http://www.jneurosci.org/content/35/19/7403.abstract
Disclaimer: These scripts are unorganized, incomplete and not as good documented as they could/should be

I put them online anyway because 
  a) I'm interested in all mistakes I made! Please tell me so I can act accordingly!
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
  1. run bs_limo_Generate(subjToRun, flagIdx, runlocal)
      This does several key things:
          a) LIMO = bs_limo_defineLIMO('subject',subjList{subj},'dataDir',data_merge);
              % Defines the LIMO Structure with all paths, options etc.

          b) LIMO = bs_limo_designMat(LIMO,'full','yes','stim',opt.stim,'interaction','full');
              % calculates the designmatrix, I added custom descriptions to all my predictors just to know which ones are actualy which
          c) limo_eeg(4) %run the actual GLM fit
          


  2. run bs_eeg_processMainGroup
      This runs the statistics for each predictor over subject (i.e. second level statistics, t-test in our case)
          a) bs_limo_generate2ndLIMO
              % Defines the second level LIMO structure
          b) bs_limo_ttestContrast, bs_limo_interactionContrast,bs_limo_interactionContrast_three_way,bs_limo_interactionContrast_four_way,
            run the actual ttests
  3. run bs_limo_tfceGrid
        This is actually very ugly, I call my custom plotting function, which internally calls the multiple comparison correction of choice (TFCE here) and saves the multiple comparison correction. As LIMO saves all permutation t-test values beforehand, the TFCE calculation is comparably quick. the plotting function is: bs_limo_display_resultsV2 interesting parts start from L137, basically I call: bs_limo_tfceCalc and then save the result
        TLDR; run bs_limo_tfceCalc
          
  4. I visualize the results using bs_limo_display_resultsV4
