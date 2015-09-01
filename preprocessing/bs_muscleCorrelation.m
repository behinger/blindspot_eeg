[testSpectra freqs] = spectopo( EEG.icaact(19,:), EEG.pnts, EEG.srate,'freqrange',[1 50]);
plot(freqs(freqs>1 &freqs<50),testSpectra(freqs>1 &freqs<50))


for k = 1:40
[spectra(k,:),freqs] = spectopo( EEG.icaact(k,:), EEG.pnts, EEG.srate,'freqrange',[1 50],'plot','off');
end
%%
for k = 1:40
w = [ones(1,28),linspace(1,0.1,72)];
R = weightedcorrs([perfectMuscle;spectra(k,(freqs>1 &freqs<50))]',w);
corrR(k) = corr(testSpectra(freqs>1 &freqs<50)',spectra(k,(freqs>1 &freqs<50))');
corrR2(k) = R(2,1);
if corrR(k) > 0.8
    figure
    title(num2str(k))
    plot(freqs(freqs>1 &freqs<50),testSpectra(freqs>1 &freqs<50))
hold on,plot(freqs(freqs>1 &freqs<50),spectra(k,freqs>1 &freqs<50))
legend('test','tobetested')
end
end
