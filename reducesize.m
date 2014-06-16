function [xr, yr] = reducesize(x, y, size)
    windowSize = floor(length(y)/(size));
    yf = filter(ones(1,windowSize)/windowSize,1,y);
    xf = filter(ones(1,windowSize)/windowSize,1,x);
    yr = yf(windowSize:windowSize:windowSize*size);
    xr = xf(windowSize:windowSize:windowSize*size);
end