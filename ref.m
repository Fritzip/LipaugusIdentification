function Imax = ref(filename,t1, t2)
    [x, fs] = audioread(filename);
    x = x(ceil(t1*fs):floor(t2*fs));
    
    % Normalization
    xmax = max(abs(x));
    x = x/xmax;

    % Signal parameters
    xlen = length(x);
    t = (0:xlen-1)/fs;

    % Define analysis and synthesis parameters
    wlen = 1024;
    h = wlen/4;
    nfft = wlen;
    
    % Bandpass filter
    x = bandpass(x,1000,6000);
    
    % Perform time-frequency analysis
    [S, ~, ~] = stft(x, wlen, h, nfft, fs);
    
    modulex = abs(S);

    im = mat2gray(modulex);
    imtotreat = im*255;
    
    m = 1;
    Imax = zeros(size(imtotreat,2),m);
    Vmax = zeros(size(imtotreat,2),m);
    for i = 1:size(imtotreat,2)
        [sortedValues, sortedIndex] = sort(imtotreat(:,i),'descend');
        maxIndex = sortedIndex(1:m);
        maxValues = sortedValues(1:m);
        Imax(i,:) = maxIndex;
        Vmax(i,:) = maxValues;
    end
    %imagesc(imtotreat)
    %set(gca,'YDir','normal')
    %hold on
    Imax(Vmax<40) = NaN;
    %plot(1:size(imtotreat,2), Imax, '*r')
end