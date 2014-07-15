function m = createteamplate(x,m)
%%%%%%%%%%%%%%%%%%%%%%%
% Initalization
%%%%%%%%%%%%%%%%%%%%%%%
lipoalgopaths;
tic
global fs nfft ovlp T F Fint

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
% Pre-treatment
%%%%%%%%%%%%%%%%%%%%%%%
tic
% Bandpass filter
xbp = bandpass(x,LOW,HIGH);

% Spectro
[S, F, T, ~] = dospectro(xbp,ovlp);

% Keep only interesting part
Sint = S(freq2co(LOW):freq2co(HIGH),:);
Fint = F(freq2co(LOW):freq2co(HIGH));

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
% Segmentation
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
%figure(1)
%plotmat(T,Fint,log(Sa)), hold on
% 
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
for i = 2:length(pks(:,1))
    %disp(i)
    
    % Get position (time) of calls
    callstart = pks(i,1);
    [tinf, tsup] = findtco(callstart, shftw, size(Safn,2)); %%%%%%%%%%% /!\ ARBITRARY CONST IN FUNCTION %%%%%%%%%%%
    
    trange = tinf:tsup;
    frange = freq2coint(LOWR):freq2coint(HIGHR);
    
    % Reduce mat to the studied call
    Sreca = Sa(frange, trange);
    Srecaf = Saf(frange, trange);

    % Find the edges (denoise)
    BWsa = findedges(Sreca);
    BWsaf = findedges(Srecaf);
    BWs = [BWsa, BWsaf];
    
    % Apply BWareaopen (matlab function)
    m = m + rmnoisepts(BWs,freq2co(220)*time2co(0.1)); %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%% taille minimale bloc en pixels 6*6 
    figure(1), plotmat(m)
end

end
