function out = co2time(varargin)
%     global fs nfft ovlp
%     out = zeros(1,length(varargin));
%     for i = 1:length(varargin)
%         out(i) = (100-ovlp)*varargin{i}*nfft/(100*fs);
%     end
    global T
    out = zeros(1,length(varargin));
    for i = 1:length(varargin)
        out(i) = T(varargin{i});
    end
end