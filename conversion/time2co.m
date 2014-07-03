function out = time2co(varargin)
%     global fs nfft ovlp
%     out = zeros(1,length(varargin));
%     for i = 1:length(varargin)
%         out(i) = round(varargin{i}*100*fs/(nfft*(100-ovlp)));
%     end
    global T
    out = zeros(1,length(varargin));
    for i = 1:length(varargin)
        ind = find(T<varargin{i});
        if isempty(ind)
            ind = 1;
        end
        out(i) = ind(end);
    end
end