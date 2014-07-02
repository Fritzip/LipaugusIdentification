function msr = computeareasum(seg, val)
%computemeasurments : Make measurments on Pi signal
%   Create a vector whose sum is equal to the area of pi signal from y =
%   val to y = max. Save the value (in msr) of the partial area at each row.
    msr = 0;
    if ~isequal(seg,0)
        indinf = find(seg.yq==val,1,'first');
        indsup = length(seg.yq);
        valsup = seg.xq(indsup);
        for i = indinf:indsup-1
            msr(end+1) = msr(end) + valsup - seg.xq(i);
        end
    end
end

