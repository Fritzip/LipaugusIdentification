%function measures = lipoalgo(x, fs_, cut, k)

%%%%%%%%%%%%%%%%%%%%%%%
%% Initalization
%%%%%%%%%%%%%%%%%%%%%%%
tic, clear, clc, close all
lipoalgopaths;

global fs nfft ovlp T F Fint Tint

% Read .wav file
cut = 2600;
k = 1;
[x, fs] = audioread('120119_071_mono3.wav',[round(cut*44100) round((cut+200)*44100)]);
%fs = fs_;

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
Sint = S(freq2co(LOW):freq2co(HIGH),:);
Fint = F(freq2co(LOW):freq2co(HIGH));
Tint = T + cut;

% Denoise signal
%Sdn = denoise(S,-20);
%Sdna = abs(Sdn);

% Prepare image
Sa = abs(Sint);
im = mat2gray(Sa)*255;

% Fuzzy Filter
h=fspecial('average');
Saf =filter2(h,Sa);

% Neighborhood and Block Processing
f = @(x) min(x(:));

%%%%%%%%%%% /!\ CONSTANTES ARBITRAIRES MAIS DÉPENDANTES DES C.I. %%%%%%%%%%%
%%%%%%%%%%% (en pixel ~ carré de 2*20)
Safn = nlfilter(Saf,[freq2co(850) time2co(0.04)],f); 

% % Plot Spectro 
% plotmat(T,Fint,log(Sa));
% % Plot Fuzzy Spectro
% plotmat(T,Fint,log(Saf));
% % Plot Fuzzy Neigborhooded Spectro
% plotmat(T,Fint,log(Safn));

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
    Srecafn = Safn(freq2coint(LOWR):freq2coint(HIGHR),i:i+shftw);
    sumit(i/pas) = sum(Srecafn(:));
end
sumit = sumit/max(sumit);

% Peak detection
delta = 0.005;         %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%
% delta correspond à la sensibilité de la détection. Plus delta diminu,
% plus on obtient de pks (peaks), et inversement. On peut alors ce servir
% de ce paramètre pour être plus ou moins sensible. 
[pks, ~] = peakdet(sumit, delta, 1:length(sumit));

%figure(1), plot(sumit), hold on, plot(pks(:,1),pks(:,2),'*r'), hold off

% Plot spectro and rectangles
figure(k)
plotmat(Tint,Fint,log(Sa)), hold on

% for i = 1:size(pks,1)
%     dec = 30; % en Hz
%     if isequal(mod(i,2),0)
%         color = 'b';
%         dec = -dec;
%     else
%         color = 'k';
%     end
%     rectangle('Position',[co2time(pks(i,1))-0.15,1100+dec,1.6,4500],...
%             'Curvature',[0.4,0.8],...
%             'LineWidth',2,...
%             'LineStyle','--',...
%             'EdgeColor',color)
%     hold on
%     text(co2time(pks(i,1))+0.6,5800,num2str(i),...
% 	'VerticalAlignment','middle',...
% 	'HorizontalAlignment','center',...
% 	'FontSize',14)
% end
% hold off
toc

%%%%%%%%%%%%%%%%%%%%%%%
% Treatment 
%%%%%%%%%%%%%%%%%%%%%%%

seg = cell(length(pks(:,1)),1);
measures = [];
ni = 1;
truth = {{'a' 'b' 'c' 'd' 'd' 'c' 'a' 'c' 'd' 'e' 'f' 'd' 'g' 'b' 'h' 'g' 'b'}... 
{'h' 'g' 'd' 'c' 'c' 'b' 'd' 'g' 'b' 'c' 'g' 'b' 'h' 'd' 'b' 'd' 'e'}...
{'c' 'Y' 'b' 'd' 'e' 'd' 'b' 'd' 'b' 'd' 'b' 'g' 'f' 'd' 'b' 'e' 'f' 'd' 'g' 'd' 'e' 'b' 'j' 'b' 'g'}...
{'a' 'c' 'j' 'b' 'a' 'g' 'j' 'b' 'e' 'd' 'j' 'b' 'g' 'a' 'j' 'b' 'e' 'b' 'j' 'a' 'b' 'j' 'g' 'a'}...
{'g' 'j' 'a' 'b' 'a' 'g' 'j' 'g' 'd' 'b' 'j' 'c' 'g' 'a' 'b' 'e' 'X' 'j' 'g' 'a' 'c' 'e' 'j' 'd' 'b'}};

template = double(rgb2gray(imread('thetemplate4.png')));
template = flipdim(template,1);
template = imresize(template,[101,78]); %%%%% CONST %%%%%%

for i = 1:length(pks(:,1))
    disp(i)
    
    % Get position (time) of calls
    [tinf, tsup] = findtco(pks(i,1), shftw, size(Safn,2)); %%%%%%%%%%% /!\ ARBITRARY CONST IN FUNCTION %%%%%%%%%%%
    
    trange = tinf:tsup;
    frange = freq2coint(LOWR):freq2coint(HIGHR);
    
    % Reduce mat to the studied call
    Sreca = Sa(frange, trange);
    Srecaf = Saf(frange, trange);
    
    new = template.*(Sreca);%+abs(min(log(Sreca(:)))));
    %new = new+abs(min(new(:)));
    %newdn = (new>0.05*max(new(:))).*new; %%%%% CONST %%%%%%

    
    % Find the edges (denoise)
    BWsa = findedges(Sreca);
    BWsaf = findedges(Srecaf);
    BWs = [BWsa, BWsaf];
    
    % Apply BWareaopen (matlab function)
    m = rmnoisepts(BWs,freq2co(220)*time2co(0.1)); %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%% taille minimale bloc en pixels 6*6
    
    figure(5)
    subplot(131), plotmat(m)
    m = m.*template;
    m(m<0.1*max(m(:))) = 0;
    subplot(132), plotmat(m)

    %%%%%%%%%%%%%%%%%%%%%%%
    % Measurments
    %%%%%%%%%%%%%%%%%%%%%%%

    % Edges detector - get pieces of curve from matrix
    [seg{i}, mout] = mat2pocs(m);
    [I, J] = ind2sub(size(mout),find(mout>0));
    
    subplot(133), plotmat(mout)
    hold on
    [x1,y1,w1] = mat3vec(mout(:,1:round((35/78)*size(Sreca,2))));
    [x2,y2,w2] = mat3vec(mout(:,round((36/78)*size(Sreca,2)):size(Sreca,2)));
    [smoothfit1, ~] = createFitSmooth2(x1, y1, w1, 0.0045); %%%%% CONST %%%%%%
    [smoothfit2, ~] = createFitSmooth2(x2, y2, w2, 0.0045); %%%%% CONST %%%%%%
    
    plot(6:35, smoothfit1(6:35).*((smoothfit1(6:35)-20)>0), 'r', 'LineWidth',2)
    plot(39:75, smoothfit2(3:75-39+3).*((smoothfit2(3:75-39+3)-20)>0), 'r', 'LineWidth',2)
    hold off
        
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
        
        % Linear Fitting
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
    value = freq2coint(1600); %%%%%%% CONST %%%%%%%
    piseq = getpisignal(seg{i},value);
    
    areasum = computeareasum(piseq,value);
    %areasumpi = computeareasum2(smoothfit1(1:0.1:30));
    %areasumhau = computeareasum2(smoothfit2(1:0.1:40));
    
    if length(areasum) > 6 && max(areasum) > 400 %%%%%%%% CONST %%%%%%%%
        
        figure(k)
        xlim([co2timeint(tinf)-1 co2timeint(tsup)+1])
        
        dec = 30; % en Hz
        if isequal(mod(i,2),0)
            color = 'b';
            dec = -dec;
        else
            color = 'k';
        end
        rectangle('Position',[co2timeint(pks(i,1))-0.15,1100+dec,1.6,4500],...
                'Curvature',[0.4,0.8],...
                'LineWidth',2,...
                'LineStyle','--',...
                'EdgeColor',color)
        hold on
        text(co2timeint(pks(i,1))+0.6,5800,truth{k}{ni},...
        'VerticalAlignment','middle',...
        'HorizontalAlignment','center',...
        'FontSize',14)

        au = (smoothfit1(1:0.1:30)-20).*((smoothfit1(1:0.1:30)-20)>0);
        ie = (smoothfit2(1:0.1:40)-20).*((smoothfit2(1:0.1:40)-20)>0);
        
        fitresults = createFitFourier2(areasum);
        measures = [measures; tinf double(uint8(truth{k}{ni})) fitresults.a0 fitresults.a1 fitresults.b1...
                    fitresults.a2 fitresults.b2 fitresults.w...
                    max(areasum) length(areasum) sum(Sreca(:)) increasesize(areasum,size(m,1)-value) au' ie'];


        
%         measures = [measures; tinf double(uint8(truth{k}{ni})) ];
            
            %max(areasum) length(areasum) sum(Sreca(:))];  
                
                
        ni = ni+1;

%     else
%         measures = [measures; i 0 0 0 0 0 0 max(areasum) length(areasum) sum(Sreca(:)) increasesize(areasum,size(m,1)-value)];
%         pks(i,:) = []; % on supprime la ligne d'un potentiel lipaugus non exploitable

    end
    
    hold off 
    
    PLOT = 1;
    if PLOT
        % Focus on the current signal
        figure(1)
        xlim([co2timeint(tinf)-1 co2timeint(tsup)+1])

        % Plot
        figure(2)
        subplot(131), plotmat(log(Sreca)); title('Raw') %Tint(trange),Fint(frange),
        hold on
        plot(6:35, smoothfit1(6:35).*((smoothfit1(6:35)-20)>0), 'k', 'LineWidth',2)
        plot(39:75, smoothfit2(3:75-39+3).*((smoothfit2(3:75-39+3)-20)>0), 'k', 'LineWidth',2)
        hold off
        subplot(132), plotmat(Tint(trange),Fint(frange), m); title('The matrix to treat') % plotmat(new)
        subplot(133), plotseg(seg{i},0,1,0,0) % plotmat(newdn)

        figure(3)
        plot(areasum,'-r'), hold on, %xlim([30 52]), ylim([400 1000])
            
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

%end

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


