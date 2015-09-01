function [EEG,muscleRej,spectRej,raw] = be_ICA_detectMuscle(EEG,p)
%% [EEG,muscleRej,spectRej] = be_ICA_detectMuscle(EEG,p)
% muscleRej is a correlation to a stereotypical muscle spectrum
% spectRej is the highfreqs/lowfreq > 2 of Jose Ossandon

perfectMuscle = [-23.91,-25.43,-26.72,-27.50,-27.92,-28.11,-28.21,-28.20,-28.02,-27.64,-27.08,-26.49,-25.96,-25.40,-24.72,-23.95,-23.25,-22.73,-22.44,-22.35,-22.31,-22.21,-22.08,-21.99,-21.89,-21.78,-21.68,-21.60,-21.56,-21.52,-21.46,-21.40,-21.32,-21.20,-21.02,-20.89,-20.83,-20.77,-20.69,-20.66,-20.67,-20.64,-20.56,-20.43,-20.36,-20.37,-20.41,-20.43,-20.43,-20.44,-20.49,-20.55,-20.57,-20.58,-20.56,-20.54,-20.54,-20.55,-20.56,-20.58,-20.59,-20.58,-20.58,-20.61,-20.64,-20.63,-20.61,-20.61,-20.59,-20.55,-20.56,-20.57,-20.54,-20.53,-20.56,-20.59,-20.59,-20.60,-20.65,-20.70,-20.70,-20.66,-20.68,-20.74,-20.77,-20.77,-20.74,-20.71,-20.68,-20.67,-20.67,-20.68,-20.70,-20.71,-20.71,-20.69,-20.69,-20.71,-20.74,-20.73];

fprintf('calculating Spectra \n')

%     [~,freqs,spectra] = spectopo( EEG.data, EEG.pnts, EEG.srate,'freqrange',[1 100],'plot','off','weights',EEG.icaweights);
ft_progress('init','etf')
for k = 1:size(EEG.icawinv,2)
%     [spectra(k,:),freqs] = spectopo( EEG.icaact(k,:), EEG.pnts,EEG.srate,'freqrange',[1 100],'plot','off','nfft',1024,'winsize',1024,'mapnorm',EEG.icawinv(k,:));
    ft_progress(k/size(EEG.icawinv,2),'Calculating Spectrum %d from %d',k,size(EEG.icawinv,2))
    [spec f] = pwelch(EEG.icaact(k,:),1024,0,1024,500);
%     spectra(k,:) = 10.^spectra(k,:);
    % w = sqrt(1./freqs);
    % weight the beginning of the sqrt-shape curve higher than the tail
%     w = [ones(1,28),linspace(1,0.1,72)];
%     R = weightedcorrs([perfectMuscle;spectra(k,(freqs>1 &freqs<50))]',w);
%     corrR(k) = R(2,1);
%     corrR(k) = corr(perfectMuscle(1:EEG.srate/500:end)',spectra(k,(freqs>1 &freqs<50))');

    specSel = spec((f>1 &f<50));
    try
        corrR(k) = corr(perfectMuscle',log10(specSel));
    
    catch
        %this was here before
        %corrR(k) = corr(perfectMuscle(1:EEG.srate/500:end)',log10(specSel));
        error('modified this file for popout2 sampling rate of 1000, why did this not work with the other project?')
    end

%     lowFreq(k) = mean(spectra(k,(freqs>1 & freqs <20)));
%     highFreq(k) = mean(spectra(k,(freqs>20 & freqs <100)));
   lowFreq(k) = median(spec(f>1&f<20));
   highFreq(k)=  median(spec(f>20&f<100));
    
end
ft_progress('close')
  badSpect = highFreq./lowFreq > 2;
raw.spect = highFreq./lowFreq;
raw.corr = corrR;
fprintf('Following Components have been marked :\n')
fprintf('%i ',find(corrR>0.7 | badSpect))
fprintf('\n Total of: %i \n',sum(corrR>0.7 | badSpect))
if ~isfield(EEG.reject,'gcompreject') || isempty(EEG.reject.gcompreject)
    EEG.reject.gcompreject = zeros(1,size(EEG.icawinv,2));
end
muscleRej = corrR>0.7;
spectRej = badSpect;
EEG.reject.gcompreject =  EEG.reject.gcompreject | badSpect | muscleRej;
end