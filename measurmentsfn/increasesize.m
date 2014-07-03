function [ vec ] = increasesize(vec,size)
last_val = vec(end);
vec = [vec repmat(last_val,1,size-length(vec))];
end

