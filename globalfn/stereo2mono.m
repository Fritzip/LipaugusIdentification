function wout = stereo2mono(win)
    if size(win,2) >= 2
        wout = (win(:,1)+win(:,2))/2;
    else
        wout = win;
    end
end