function ffft = getffft(wffta)
    global fs
    ffft = (0:length(wffta)-1)*(fs/(2*(length(wffta)-1)));
end