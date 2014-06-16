function plotoscillo(t, w)
    global fs nbech
    
    plot(t, w)
    %yPos = std(w)*2;
    %hold on
    %plot(get(gca,'xlim'), [yPos yPos], 'r');
    %hold off
    xlabel('Time (s)')
    ylabel('Amplitude')
    title('Oscillo')
    xlim([0 nbech/fs])
end