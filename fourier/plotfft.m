function plotfft(f, wffta)
    plot(f, wffta)
    xlabel('Frequency (Hz)')
    ylabel('Amplitude')
    title('Fourier Transform') 
    xlim([0 10000])
end