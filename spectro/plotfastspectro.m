function plotfastspectro(w)
    global fs nfft
    specgram(w,nfft,fs)
end