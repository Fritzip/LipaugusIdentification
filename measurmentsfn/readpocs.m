function [x, y] = readpocs(pocs)
    data = sortrows(pocs,1);
    x = data(:,1);
    y = data(:,2);
end