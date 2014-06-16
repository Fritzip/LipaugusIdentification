function plotmat(varargin)
    if isequal(size(varargin,2),1)
        imagesc(varargin{1});
    elseif isequal(size(varargin,2),3)
        imagesc(varargin{1}, varargin{2}, varargin{3});
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
    end
    set(gca,'YDir','normal');
    %title('Spectro')
end