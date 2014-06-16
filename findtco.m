function [tinf, tsup] = findtco(val, shftw)
    tinf = val-time2co(0.15);
    tsup = val+shftw;
    if tinf < 1
        tinf = 1;
    end
end