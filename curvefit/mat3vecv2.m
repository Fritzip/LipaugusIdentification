function [X, Y, W] = mat3vecv2(mat)
% Transform mat to 3 vectors x,y (positions) and w the weighted vector
% -------------------------------------------------------------------------
[Y, X] = ind2sub(size(mat),find(mat>0));
W = mat(mat>0);
end