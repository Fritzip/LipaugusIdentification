function [x, y, w] = mat3vec(m)
% Transform m to 3 vectors x,y (positions) and w the weighted vector
% -------------------------------------------------------------------------
    n = size(m,1);
    y = zeros(size(m,2),n);
    w = zeros(size(m,2),n);

    for j = 1:size(m,2)
        [sortedValues, sortedIndex] = sort(m(:,j),'descend');
        y(j,:) = sortedIndex(1:n);
        w(j,:) = sortedValues(1:n);
    end

    y = y'; y = y(:);
    w = w'; w = w(:);
    x = repmat((1:size(m,2))',1,n)'; x = x(:)';
end