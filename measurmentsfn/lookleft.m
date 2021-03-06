function [xnew, ynew, bool] = lookleft(x, y, m)
% Called in case of the following situation to find B px when you are in A.
%      ___
%     |i,j|
%     |_B_| ___
%          |x,y|
%          |_A_|
%   
    xnew = x;
    ynew = y;
    zerotoone = 0;
    onetozero = 0;
    dec = 0;
    while (~zerotoone || ~onetozero) && dec <= 3  %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%
        [xnew, ynew] = checkco(xnew - 1, ynew, m);
        if ~onetozero && ~zerotoone && ~isequal(m(ynew,xnew),0)
            zerotoone = 1;
        end
        if ~onetozero && zerotoone && isequal(m(ynew,xnew),0)
            onetozero = 1;
        end
        dec = dec + 1;
    end
    if zerotoone && onetozero
        bool = 1;
        xnew = xnew + 1;
    else
        bool = 0;
    end
end