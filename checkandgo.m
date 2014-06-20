function [xnew, ynew, bool] = checkandgo(x, y, m, dir)
    [xnew, ynew] = checkco(x + dir, y, m);
%     if ~isequal(xnew,x) || ~isequal(ynew,y)
%         % on a trouver un bord
%         bool = 1;
%         return
    if isequal(m(ynew,xnew),0) && isequal(dir,-1)
        bool = 1;
        xnew = x;
        ynew = y;
    elseif ~isequal(m(ynew,xnew),0) && isequal(dir,1)
        bool = 1;
    else
        bool = 0;
    end
end
