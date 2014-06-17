%%%%%%%%%%%%%%%%%%%%%%%
%% Initalization
%%%%%%%%%%%%%%%%%%%%%%%
tic, clear, clc, close all
lipoalgopaths;

global fs nfft ovlp T F

% Read .wav file
[x, fs] = audioread('120119_071_mono3.WAV',[round(2700*44100) round(2900*44100)]);

% Stereo to mono
x = stereo2mono(x);

% Normalization
xmax = max(abs(x));
x = x/xmax;

% Signal parameters
xlen = length(x);
t = (0:xlen-1)/fs;

% Define analysis and synthesis parameters
wlen = 1024;
ovlp = 25;
hop = round(ovlp*wlen/100);
nfft = wlen;

% Constants
HIGH = 6000;
LOW = 1000;
HIGHR = 5600;
LOWR = 1300;

toc
%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-treatment
%%%%%%%%%%%%%%%%%%%%%%%
tic
% Bandpass filter
xbp = bandpass(x,LOW,HIGH);

% Spectro
[S, F, T, ~] = dospectro(xbp,ovlp);

% Keep only interesting part
S = S(freq2co(LOW):freq2co(HIGH),:);
F = F(freq2co(LOW):freq2co(HIGH));

% Denoise signal
%Sdn = denoise(S,-20);
%Sdna = abs(Sdn);

% Prepare image
Sa = abs(S);
im = mat2gray(Sa)*255;

% Fuzzy Filter
h=fspecial('average');
Saf =filter2(h,Sa);

% Neighborhood and Block Processing
f = @(x) min(x(:));
Safn = nlfilter(Saf,[20 2],f); %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%

% % Plot Spectro 
% plotmat(T,F,log(Sa));
% % Plot Fuzzy Spectro
% plotmat(T,F,log(Saf));
% % Plot Fuzzy Neigborhooded Spectro
% plotmat(T,F,log(Safn));

toc
%%%%%%%%%%%%%%%%%%%%%%%
%% Segmentation
%%%%%%%%%%%%%%%%%%%%%%%

tic
% Sum of all intensity
shftw = time2co(1.2); % Pie + Hau average lengths %%%%%% /!\ CONST %%%%%%

pas = 1;
sumit = zeros(ceil((size(Safn,2)-shftw)/pas),1);
for i = pas:pas:(size(Safn,2)-shftw)
    sumit(i/pas) = sum(sum(abs(Safn(freq2co(LOWR):freq2co(HIGHR),i:i+shftw))));
end

% Peak detection
delta = 20;         %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%
[maxtab, mintab]=peakdet(sumit, delta, 1:length(sumit));

% Plot spectro and rectangles
figure(1)
plotmat(T,F,log(Sa)), hold on

for i = 1:size(maxtab,1)
    dec = 30;
    if isequal(mod(i,2),0)
        color = 'b';
        dec = -dec;
    else
        color = 'k';
    end
    rectangle('Position',[co2time(maxtab(i,1))-0.2,1100+dec,1.6,4500],...
            'Curvature',[0.4,0.8],...
            'LineWidth',2,...
            'LineStyle','--',...
            'EdgeColor',color)
    hold on
end
hold off
toc
%%%%%%%%%%%%%%%%%%%%%%%
%% Treatment 
%%%%%%%%%%%%%%%%%%%%%%%

i = 31;
%for i = 2:length(maxtab(:,1))
    % Get position (time) of calls
    val = maxtab(i,1);
    [tinf, tsup] = findtco(val, shftw); %%%%%%%%%%% /!\ ARBITRARY CONST IN FUNCTION %%%%%%%%%%%
    trange = tinf:tsup;
    frange = freq2co(LOWR):freq2co(HIGHR);
    
    % Reduce mat to the studied call
    Sreca = Sa(frange, trange);
    Srecaf = Saf(frange, trange);

    % Find the edges (denoise)
    m = findedges(Sreca, Srecaf); %%%%%%%%%%% /!\ ARBITRARY CONST IN FUNCTION %%%%%%%%%%%
    
    % Transform m to 3 vectors x,y (positions) and w the weighted vector
    n = size(m,1);
    y = zeros(size(m,2),n);
    w = zeros(size(m,2),n);

    for j = 1:size(m,2)
        [sortedValues, sortedIndex] = sort(m(:,j),'descend');
        y(j,:) = sortedIndex(1:n);
        w(j,:) = sortedValues(1:n);
    end

    y = y'; y = y(:);
    w = w'; w = w(:);
    x = repmat(T(trange)',1,n)'; x = x(:)';

    % Plot
    figure(2)
    subplot(131), plotmat(T(trange),F(frange),log(Sreca)); title('Raw')
    subplot(132), plotmat(T(trange),F(frange),m); title('The matrix to treat')
    subplot(133), plotmat(T(trange),F(frange),m); title('What is expected') % plot(x, y, '*b')
    
%     waitforbuttonpress
% end


%%%%%%%%%%%%%%%%%%%%%%%
%% Measurments
%%%%%%%%%%%%%%%%%%%%%%%
tic
for k = length(maxtab(:,1))
    Srec = Safn(freq2co(LOWR):freq2co(HIGHR),maxtab(k,1)-time2co(0.1):maxtab(k,1)+shftw);
    for i = 1:size(Srec,1)
        [pks, locs, over] = getpeakseparated(diff(Srec(i,:)),time2co(0.05),5);
        pks(~logical(over)) = [];
        locs(~logical(over)) = [];

        subplot(211), plotmat(log(Srec));
        subplot(212), plot(diff(Srec(i,:)))
        hold on
        plot(locs,pks, '*r')
        hold off
        waitforbuttonpress;
    end
    for j = 1:size(Srec,2)
        
    end
    break
end
toc


%% getpeakseparated
[pks, locs, over] = getpeakseparated(sumit,1,0,1.3*min(sumit));
pks(~logical(over)) = [];
locs(~logical(over)) = [];

figure
plot(sumit)
hold on
plot(locs,pks, '*r')


%%

BW = edge(C,'zerocross');
plotmat(T,F,BW);

%%
[pks, locs, over] = getpeakseparated(Sa(freq2co(4000),:),15,0,0.5);
pks(~logical(over)) = [];
locs(~logical(over)) = [];

plot(locs, pks, '*r');

% Find m max values by row
m = 10;
Imax = zeros(size(im,2),m);
Vmax = zeros(size(im,2),m);
for i = 1:size(im,2)
    [sortedValues, sortedIndex] = sort(im(:,i),'descend');
    maxIndex = sortedIndex(1:m);
    maxValues = sortedValues(1:m);
    Imax(i,:) = maxIndex;
    Vmax(i,:) = maxValues;
end
%Imax(Vmax<5) = NaN;
%plot(T, Imax*fs/nfft, '*b')

%%

%%
[C, I] = max(im);
I(C<45) = NaN;
I = fixgaps(I);
I(isnan(I)) = [];

plot(I)

lenI = length(I);

% Assume the sampling interval of the signal Z is 1/32 seconds. Construct
% an appropriate time axis and plot the signal
time = linspace(0,h*lenI/fs,lenI);
figure;
plot(time,I); grid on;


% From the figure, we see that one pattern is centered around 20
% seconds and the second pattern is centered around 48 seconds.
% The pattern centered around 20 seconds has a duration of 8 seconds while
% the second pattern has a duration of 4 seconds.

% We analyze the signal by computing the CWT coefficients of Z using the
% admissible wavelet we constructed to approximate the basic form F. 

% stepSIG = 1/fs;
% stepWAV = 1/fs;
% wname = 'lipo';
% scales  = (1:2*lenI)*stepSIG;
% WAV = {wname,stepWAV};
% SIG = {I,stepSIG};
% figure;
% cwt(SIG,scales,WAV,'scalCNT'); grid


%%
I = I/max(I);
I1 = I;%(3500:4500); 
scales = 1:256;
wname = 'lipo';

% Continuous wavelet transform (CWT).
figure;
c = cwt(I1,scales,wname,'scalCNT'); grid


max(abs(c(:)))



% Detect the maximum of the absolute value of coefficients.
lenSIG = length(I);
%positions = (0:lenSIG-1)*stepSIG;    
%fprintf('Instant 1:  %4.2f\n',positions(round(middle)))
%fprintf('Instant 2:  %4.2f\n',positions(round(middle2)))



%% Max of each column
[C, I] = max(imtotreat);

for i = 3:size(imtotreat,2)-3
    
    leftneighbor = [i-2 I(i-2); i-1 I(i-1)];
    rightneighbor = [i+1 I(i+1); i+2 I(i+2)];
    distal = mean(pdist2(leftneighbor, [i I(i)], 'euclidean'));
    distar = mean(pdist2(rightneighbor, [i I(i)], 'euclidean'));
    
%     distw = zeros(1,3);
%     for j = 1:3
%         distw(j) = pdist([neighbor(j,:); neighbor(j+1,:)],'euclidean');
%     end
%     distw = mean(distw);
    
    distwl = pdist([leftneighbor(1,:); leftneighbor(2,:)],'euclidean');
    distwr = pdist([rightneighbor(1,:); rightneighbor(2,:)],'euclidean');
    
    if (distal-distwl)*100/distwl > 50 && (distar-distwr)*100/distwr > 50
        I(i) = I(i-1);
    end
end

% 
% meandist = pdist([1 I(1); 2 I(2)],'euclidean');
% alldist = [];
% SAUT = 0;
% for i = 3:length(I)-3
%     SAUT = 0;
%     currdist = pdist([i-1 I(i-1); i I(i)],'euclidean');
%     if (currdist-meandist)*100/meandist > 25 && i>20
%         
%     else
%         alldist = [alldist currdist];
%         meandist = mean(alldist);
%     end
% end

%% Only high amplitude
I(C<35) = NaN;
imagesc(imtotreat)
set(gca,'YDir','normal')
hold on
plot(I,'-*r')
%plot(maxes(:,1),maxes(:,2), '*r')
%imshow(imtotreat)
%figure
%plot(diff(I))
%noise = I(1:50);
%stdnoise = std(noise)
%%
plot(I, '*g')
%%
% resynthesis of the original signal
stftb = denoise(stft,0);
[x_istft, t_istft] = istft(stftb, h, nfft, fs);

surf(log(abs(stft(:,200:1500))),'edgecolor','none'); axis tight; view(0,90);
figure
surf(log(abs(stftb(:,200:1500))),'edgecolor','none'); axis tight; view(0,90);
% plot the original signal
%figure(1)
%subplot(211), plot(t(1:1000), x(1:1000), 'b')
%grid on
%xlim([0 max(t)])
%ylim([-1.1 1.1])
%set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
%xlabel('Time, s')
%ylabel('Normalized amplitude')
%title('Original and reconstructed signal')

% plot the resynthesized signal 
%hold on
%subplot(212), plot(t_istft(1:1000), x_istft(1:1000), 'r')
%legend('Original signal', 'Reconstructed signal')
%xlim([0 max(t)])
%ylim([-1.1 1.1])

%% Wavelet
% Ref
F = ref('A.WAV',184.95,186.425);
F([27 30 94 95 129:136]) = NaN;
F = F(15:212);
F = fixgaps(F);
F = F/max(F);
F = F';
X = linspace(-1,1,length(F));
[fitresult, gof] = createFit(X, F);
%plot(X,F)
%hold on
%plot(X, fitresult(X))

fitfit = fitresult(X)';
[Y,X] = pat2cwav(fitfit,'orthconst',0,2);
figure;
plot(X,fitfit/max(abs(fitfit)),'b'); 
hold on; 
plot(X,Y/max(abs(Y)),'r');
title('Form to detect (b) and adapted Wavelet (r)')

% % Save the adapted wavelet and add it to the toolbox 
locdir = cd;
cd(tempdir);
save adp_FRM1 X Y
wavemngr('add','Lipogus','lipo',4,'','adp_FRM1.mat',[0 1]);
addpath(tempdir,'-begin');
cd(locdir);

