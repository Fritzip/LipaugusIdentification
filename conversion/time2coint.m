function out = time2coint(varargin)
    global Tint
    out = zeros(1,length(varargin));
    for i = 1:length(varargin)
        ind = find(Tint<varargin{i});
        if isempty(ind)
            ind = 1;
        end
        out(i) = ind(end);
    end
end