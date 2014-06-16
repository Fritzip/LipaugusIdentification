function out = freq2co(varargin)
%     global fs nfft
%     out = zeros(1,length(varargin));
%     for i = 1:length(varargin)
%         out(i) = round(varargin{i}*nfft/fs);
%     end

    global F
    out = zeros(1,length(varargin));
    for i = 1:length(varargin)
        ind = find(F<varargin{i});
        out(i) = ind(end);
    end
end