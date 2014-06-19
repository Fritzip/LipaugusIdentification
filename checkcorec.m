function [x, y, w, h] = checkcorec(x, y, w, h, m)
% Given coordinates (defining a rectangle area) and a matrix, CHECKCO will
% adapt x, y, w, h to define corrects coordinates of m.
% -------------------------------------------------------------------------
    [mx, my] = size(m);
    if x < 1
        x = 1;
    end
    if y < 1
        y = 1;
    end
    if x + w > my
        w = my - x;
    end
    if y + h > mx
        h = mx - y;
    end
end