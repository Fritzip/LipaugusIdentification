function [seg, tabout] = permutesegitem(seg, value)
    tab = [];
    tabout = 0;
    if ~isempty(seg)
        for i = 1:length(seg)
            ind = find(seg{i}.yq==value, 1, 'first');
            if isempty(ind)
                x = 0;
            else
                x = seg{i}.xq(ind);
            end
            tab = [tab; i x];
        end
        tab = sortrows(tab,2);
        rowi = find(tab(:,2)==0,1,'last');
        if ~isempty(rowi)
            tab = tab([rowi+1:size(tab,1) 1:rowi],:);
        end
        tabout = tab(1,2);
        seg = seg(:,tab(:,1)');
    end
end