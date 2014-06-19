function [x, y] = checkco(x, y, m)
% Given coordinates and a matrix, CHECKCO will
% adapt x, y, to define corrects coordinates of m.
% -------------------------------------------------------------------------
    [mx, my] = size(m);
    if x < 1
        x = 1;
    end
    if y < 1
        y = 1;
    end
    if x > my
        x = my;
    end
    if y > mx
        y = mx;
    end
end