function Hd = designfilter(low, high)
    global fs
    
    Fstop1 = low-300;   % First Stopband Frequency
    Fpass1 = low;   % First Passband Frequency
    Fpass2 = high;   % Second Passband Frequency
    Fstop2 = high+300;   % Second Stopband Frequency
    Astop1 = 60;     % First Stopband Attenuation (dB)
    Apass  = 1;      % Passband Ripple (dB)
    Astop2 = 60;     % Second Stopband Attenuation (dB)
    Fs     = fs;  % Sampling Frequency

    h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1, ...
                         Fpass2, Fstop2, Astop1, Apass, Astop2, Fs);

    Hd = design(h, 'equiripple');
end