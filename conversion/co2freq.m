function out = co2freq(varargin)
%     global fs nfft
%     out = zeros(1,length(varargin));
%     for i = 1:length(varargin)
%         out(i) = varargin{i}*fs/nfft;
%     end
    global F
    out = zeros(1,length(varargin));
    for i = 1:length(varargin)
        out(i) = F(varargin{i});
    end
end