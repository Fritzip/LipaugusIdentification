function [xnew, ynew, bool] = hcheckandgo(x, y, m, dir)
    [xnew, ynew] = checkco(x + dir, y, m);
    if ~isequal(xnew,x+dir)
        % x edge of the matrix
        bool = 2;
    elseif isequal(m(ynew,xnew),0) && isequal(dir,-1)
        % found a junction while going left
        bool = 1;
        xnew = x;
        ynew = y;
    elseif ~isequal(m(ynew,xnew),0) && isequal(dir,1)
        % found a junction while going right
        bool = 1;
    else
        % no junction found
        bool = 0;
    end
end
