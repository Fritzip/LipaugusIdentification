function [tinf, tsup] = findtco(val, shftw, my)
% Given a coordinate (val) and a size of shifting windows (shftw), it
% returns begin and end coordinates.
%--------------------------------------------------------------------------
    tinf = val-time2co(0.15); % slight shift to the left
    tsup = val+shftw;
    if tinf < 1
        tinf = 1;
        tsup = tinf+shftw+time2co(0.15);
    end
    if tsup > my
        tsup = my;
        tinf = tsup-shftw-time2co(0.15);
    end
end