function [x, y] = checkco(x, y, m)
% Given coordinates and a matrix, CHECKCO will
% adapt x, y, to define corrects coordinates of m.
% -------------------------------------------------------------------------
    [mx, my] = size(m);
    for i = 1:length(x)
        if x(i) < 1
            x(i) = 1;
        end
        if y(i) < 1
            y(i) = 1;
        end
        if x(i) > my
            x(i) = my;
        end
        if y(i) > mx
            y(i) = mx;
        end
    end
end