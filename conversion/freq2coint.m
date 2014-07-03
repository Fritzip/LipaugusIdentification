function out = freq2coint(varargin)
    global Fint
    out = zeros(1,length(varargin));
    for i = 1:length(varargin)
        ind = find(Fint<varargin{i});
        if isempty(ind)
            ind = 1;
        end
        out(i) = ind(end);
    end
end