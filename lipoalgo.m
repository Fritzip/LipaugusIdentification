%%%%%%%%%%%%%%%%%%%%%%%
%% Initalization
%%%%%%%%%%%%%%%%%%%%%%%
tic, clear, clc, close all
lipoalgopaths;

global fs nfft ovlp T F

% Read .wav file
%[x, fs] = audioread('lipo.WAV');
[x, fs] = audioread('120119_071_mono3.wav',[round(2600*44100) round(2630*44100)]);

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

% Constants (Hz)
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

%%%%%%%%%%% /!\ CONSTANTES ARBITRAIRES MAIS DÉPENDANTES DES C.I. %%%%%%%%%%%
%%%%%%%%%%% (en pixel ~ carré de 2*20)
Safn = nlfilter(Saf,[freq2co(1800) time2co(0.04)],f); 

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
shftw = time2co(1.2); % Pi + Hau average lengths %%%%%% CONST %%%%%%

pas = 1; % en pixel
sumit = zeros(ceil((size(Safn,2)-shftw)/pas),1);
for i = pas:pas:(size(Safn,2)-shftw)
    sumit(i/pas) = sum(sum(abs(Safn(freq2co(LOWR):freq2co(HIGHR),i:i+shftw))));
end

% Peak detection
delta = 20;         %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%
% delta correspond à la sensibilité de la détection. Plus delta diminu,
% plus on obtient de pks (peaks), et inversement. On peut alors ce servir
% de se paramètre pour être plus ou moins sensible. 
[pks, ~] = peakdet(sumit, delta, 1:length(sumit));

%figure(1), plot(sumit), hold on, plot(pks(:,1),pks(:,2),'*r'), hold off

% Plot spectro and rectangles
figure(1)
plotmat(T,F,log(Sa)), hold on

for i = 1:size(pks,1)
    dec = 30; % en Hz
    if isequal(mod(i,2),0)
        color = 'b';
        dec = -dec;
    else
        color = 'k';
    end
    rectangle('Position',[co2time(pks(i,1))-0.15,1100+dec,1.6,4500],...
            'Curvature',[0.4,0.8],...
            'LineWidth',2,...
            'LineStyle','--',...
            'EdgeColor',color)
    hold on
    text(co2time(pks(i,1))+0.6,5700,[num2str(i)],...
	'VerticalAlignment','middle',...
	'HorizontalAlignment','center',...
	'FontSize',14)
end
hold off
toc

%%%%%%%%%%%%%%%%%%%%%%%
%% Treatment 
%%%%%%%%%%%%%%%%%%%%%%%

seg = cell(length(pks(:,1)),1);
measures = [];

for i = 1:length(pks(:,1))
    disp(i)
    
    % Get position (time) of calls
    callstart = pks(i,1);
    [tinf, tsup] = findtco(callstart, shftw, size(Safn,2)); %%%%%%%%%%% /!\ ARBITRARY CONST IN FUNCTION %%%%%%%%%%%
    
    trange = tinf:tsup;
    frange = freq2co(LOWR):freq2co(HIGHR);
    
    % Reduce mat to the studied call
    Sreca = Sa(frange, trange);
    Srecaf = Saf(frange, trange);

    % Find the edges (denoise)
    BWsa = findedges(Sreca);
    BWsaf = findedges(Srecaf);
    BWs = [BWsa, BWsaf];
    
    % Apply BWareaopen (matlab function)
    m = rmnoisepts(BWs,35); %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%% taille minimale bloc en pixels

%     m1 = m(:,1:40);
%     m2 = m(:,35:55);
%     m3 = m(:,45:size(m,2));

    %%%%%%%%%%%%%%%%%%%%%%%
    % Measurments
    %%%%%%%%%%%%%%%%%%%%%%%

    % Edges detector - get pieces of curve from matrix
    seg{i} = mat2pocs(m);

    for j = 1:size(seg{i},2)
        [seg{i}{j}.xi, seg{i}{j}.yi] = readpocs(seg{i}{j}.data);

        % Normal Smoothing
%         seg{i}{j}.yy = smooth(seg{i}{j}.xi, seg{i}{j}.yi, 0.3, 'rloess'); 
%         seg{i}{j}.xs = seg{i}{j}.xi;
%         seg{i}{j}.ys = seg{i}{j}.yy;
        
        % Reverse Smoothing
        [seg{i}{j}.yy, seg{i}{j}.ind] = sort(seg{i}{j}.yi);
        seg{i}{j}.xx = smooth(seg{i}{j}.yi,seg{i}{j}.xi,0.3,'rloess'); 
        seg{i}{j}.xs = seg{i}{j}.xx(seg{i}{j}.ind);
        seg{i}{j}.ys = seg{i}{j}.yy;
        
        % Lineare Fitting
        fit = createFitLin(seg{i}{j}.xi, seg{i}{j}.yi);
        seg{i}{j}.crease = sign(fit.p1);

        % Barycentre
        sortpoc = sortrows(seg{i}{j}.data,1);
        xyb = sortpoc(round(size(seg{i}{j}.data,1)/2),:);
        seg{i}{j}.xb = xyb(1);
        seg{i}{j}.yb = xyb(2);

        % Raw Interpolation
        seg{i}{j}.yq = min(seg{i}{j}.yi):max(seg{i}{j}.yi);
        seg{i}{j}.xq = interp1(seg{i}{j}.yi, seg{i}{j}.xi, seg{i}{j}.yq);
        
        % Smooth Interpolation
%         seg{i}{j}.yq = min(seg{i}{j}.yy):max(seg{i}{j}.yy);
%         seg{i}{j}.xq = interp1(seg{i}{j}.yy, seg{i}{j}.xx(seg{i}{j}.ind),seg{i}{j}.yq);
    end
    
    % Make measurments on Pi signal
    value = 15; %%%%%%% CONST %%%%%%%
    piseq = getpisignal(seg{i},value);
    
    areasum = computeareasum(piseq,value);
    if length(areasum) > 6 && max(areasum) > 400 %%%%%%%% CONST %%%%%%%%
        fitresults = createFitFourier2(areasum);
        measures = [measures; i fitresults.a0 fitresults.a1 fitresults.b1...
                    fitresults.a2 fitresults.b2 fitresults.w...
                    max(areasum) length(areasum) sum(Sreca(:)) increasesize(areasum,size(m,1)-value)];
    end
    
    PLOT = 1;
    if PLOT
        % Center the current signal
        figure(1)
        xlim([co2time(tinf)-3 co2time(tsup)+3])

        % Plot
        figure(2)
        subplot(131), plotmat(T(trange),F(frange),log(Sreca)); title('Raw')
        subplot(132), plotmat(m); title('The matrix to treat')
        subplot(133), plotseg(seg{i},0,1,0,0)

        figure(3)
        plot(areasum,'-r'), hold on, %xlim([30 52]), ylim([400 1000])

%         figure(5),
%         subplot(131), plotmat(m1)
%         subplot(132), plotmat(m2)
%         subplot(133), plotmat(m3)
    
        % Press key to continue
        a = 1;
        while a
            a = ~waitforbuttonpress;
        end

        figure(3)
        plot(areasum,'b'), hold on, %xlim([30 52]), ylim([400 1000])
    end
end
hold off
csvwrite('measures.csv',measures)

%%
        
%     out = zeros(size(m)*100);

    %         [xx100, yy100] = checkco(round(xx(ind)*100), round(yy*100), out);
    %         out(sub2ind(size(out),yy100,xx100)) = 1;

    
%     seuil = [10 20 30 50 60 70];
% 
%     out = imresize(out, size(m));
%     out = out./max(out(:));
% 
%     delta = 0.1;         %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%
%     for j = 1:length(seuil)
%         pks = peakdet(out(seuil(j),:), delta, 1:length(out(seuil(j),:)));
%     end
%     
%     figure(4), plot(out(40,:)), hold on
%     if ~isempty(uie)
%         plot(uie(:,1),uie(:,2),'*r')
%     end
%     hold off
% 
%     figure(3), plotmat(out)



