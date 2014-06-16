function wout = bandpass(win,low,high)

    % Artisanal method    
%     bf = @(f) floor(f*nbech/fs)+1;
%     bfr = @(f) floor((fs-f)*nbech/fs);
%     
%     [wfft, ~] = dofft(win);
% 
%     wfft(bf(0):bf(low)) = 0;
%     wfft(bf(high):bf(fs/2)) = 0;
%     wfft(bfr(low):bfr(0)) = 0;
%     wfft(bfr(fs/2):bfr(high)) = 0;
%     
%     wout = ifft(wfft,'symmetric');
    
    % Filter designed properly
    Hd = designfilter(low,high);
    wout = filter(Hd, win);
end