function [S,F,T,P] = dospectro(w,ovlp)
    global fs nfft
    [S,F,T,P] = spectrogram(w,blackman(nfft),floor(nfft*ovlp/100),nfft,fs,'yaxis');
end