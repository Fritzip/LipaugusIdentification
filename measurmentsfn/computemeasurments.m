function msr = computemeasurments(seg)
%computemeasurments Summary of this function goes here
%   Detailed explanation goes here
msr = 0;
indinf = find(seg.yq==50,1,'first');
indsup = length(seg.yq);
valsup = seg.xq(indsup);
for i = indinf:indsup-1
    msr(end+1) = msr(end) + valsup - seg.xq(i);
end

end

