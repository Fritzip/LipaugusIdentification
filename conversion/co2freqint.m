function out = co2freqint(varargin)
    global Fint
    out = zeros(1,length(varargin));
    for i = 1:length(varargin)
        out(i) = Fint(varargin{i});
    end
end