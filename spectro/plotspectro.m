function plotspectro(w, ovlp)
    global fs nfft
    spectrogram(w,blackman(nfft),floor(nfft*ovlp/100),nfft,fs,'yaxis');
end