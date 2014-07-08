function out = co2timeint(varargin)
    global Tint
    out = zeros(1,length(varargin));
    for i = 1:length(varargin)
        out(i) = Tint(varargin{i});
    end
end