function [x, y, w] = readpocs(pocs)
    data = sortrows(pocs,1);
    x = data(:,1);
    y = data(:,2);
    w = 1;
end