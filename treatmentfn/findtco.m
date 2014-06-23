function [tinf, tsup] = findtco(val, shftw)
% Given a coordinate (val) and a size of shifting windows (shftw), it
% returns begin and end coordinates.
%--------------------------------------------------------------------------
    tinf = val-time2co(0.15);
    tsup = val+shftw;
    if tinf < 1
        tinf = 1;
    end
end