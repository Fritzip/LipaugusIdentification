function [wfft wffta] = dofft(w)
    wfft = fft(w);
    wffta = abs(wfft);
    if mod(length(w),2)==0
        wlen = length(w)/2;
    else
        wlen = (length(w)-1)/2;
    end
    wffta = wffta(1:wlen);  %Discard Half of Points
end