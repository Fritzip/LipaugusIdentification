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
% figure, plotmat(T,Fint,log(Saf));
% % Plot Fuzzy Neigborhooded Spectro
% figure, plotmat(T,Fint,log(Safn));

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
% plus on obtient de pks (peaks), et inversement. On peut alors se servir
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

segpi = cell(length(pks(:,1)),1);
seghau1 = cell(length(pks(:,1)),1);
seghau2 = cell(length(pks(:,1)),1);

measures = [];
ni = 1;
truth = {{'a' 'b' 'c' 'd' 'd' 'c' 'a' 'c' 'd' 'e' 'f' 'd' 'g' 'b' 'h' 'g' 'b'}... 
{'h' 'g' 'd' 'c' 'c' 'b' 'd' 'g' 'b' 'c' 'g' 'b' 'h' 'd' 'b' 'd' 'e'}...
{'c' 'Y' 'b' 'd' 'e' 'd' 'b' 'd' 'b' 'd' 'b' 'g' 'f' 'd' 'b' 'e' 'f' 'd' 'g' 'd' 'e' 'b' 'j' 'b' 'g'}...
{'a' 'c' 'j' 'b' 'a' 'g' 'j' 'b' 'e' 'd' 'j' 'b' 'g' 'a' 'j' 'b' 'e' 'b' 'j' 'a' 'b' 'j' 'g' 'a'}...
{'g' 'j' 'a' 'b' 'a' 'g' 'j' 'g' 'd' 'b' 'j' 'c' 'g' 'a' 'b' 'e' 'X' 'j' 'g' 'a' 'c' 'e' 'j' 'd' 'b'}};

templatepi = readtemplate('templatepi.png', [101 78]);
templatehau1 = readtemplate('templatehau1.png', [101 78]);
templatehau2 = readtemplate('templatehau2.png', [101 78]);

% Iterate over all segmentation results and treat the rectangle as a signal
for i = 1:length(pks(:,1))
    disp(i)
    
    % Get position (time) of calls
    [tinf, tsup] = findtco(pks(i,1), shftw, size(Safn,2)); %%%%%%%%%%% /!\ ARBITRARY CONST IN FUNCTION %%%%%%%%%%%
    
    trange = tinf:tsup;
    frange = freq2coint(LOWR):freq2coint(HIGHR);
    
    % Reduce mat to the studied call
    Sreca = Sa(frange, trange);
    Srecaf = Saf(frange, trange);
    
    %new = template.*(Sreca);%+abs(min(log(Sreca(:)))));
    %new = new+abs(min(new(:)));
    %newdn = (new>0.05*max(new(:))).*new; %%%%% CONST %%%%%%
    
    % Find the edges (denoise)
    BWsa = findedges(Sreca);
    BWsaf = findedges(Srecaf);
    BWs = [BWsa, BWsaf];
    
    % Merge each edges matrix
    m = mergeBWs(BWs);
    %m = rmnoisepts(BWs,freq2co(220)*time2co(0.1)); %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%% taille minimale bloc en pixels 6*6
    
    % Separate matrix into upper and lower part
    mpi = treatmatrix(m,templatepi); %%%%%%%%%%% /!\ ARBITRARY CONST IN FUNCTION %%%%%%%%%%%
    mhau1 = treatmatrix(m,templatehau1);
	mhau2 = treatmatrix(m,templatehau2);
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Measurments
    %%%%%%%%%%%%%%%%%%%%%%%

    % Edges detector - get pieces of curve from matrix
    [segpi{i}, mpiedg] = mat2pocs(mpi);
    [seghau1{i}, mhau1edg] = mat2pocs(mhau1);
    [seghau2{i}, mhau2edg] = mat2pocs(mhau2);
    figure(9), plotmat(mhau2edg)
    % Convert parts into vectors and weights
    [X1, Y1, W1] = mat3vec(mpiedg);
    [X2, Y2, W2] = mat3vec(mhau2edg);
    
    %%%%%%%%% En fonction !!!
    try
        [smoothfit1, ~] = createFitSmooth2(X1, Y1, W1, 0.02, []); %%%%% CONST %%%%%%

        fdata = feval(smoothfit1,X1);
        indices = abs(fdata - Y1) > 3;
        outliers = excludedata(X1,Y1,'indices',indices);
        smoothfit12 = createFitSmooth2(X1,Y1,W1,0.02,outliers);

        [smoothfit2, ~] = createFitSmooth2(X2, Y2, W2, 0.02, []); %%%%% CONST %%%%%%

        fdata = feval(smoothfit2,X2);
        indices = abs(fdata - Y2) > 3;
        outliers = excludedata(X2,Y2,'indices',indices);
        smoothfit22 = createFitSmooth2(X2,Y2,W2,0.02,outliers);
        
        pisign = sign(diff(smoothfit12(min(X1):max(X1))));
        hau2sign = sign(diff(smoothfit22(min(X2):max(X2))))*-1;
        
        inf1 = min(X1) + find(pisign==1,1,'first') - 1;
        sup1 = min(X1) + find(pisign==1,1,'last');
        inf2 = min(X2) + find(hau2sign==1,1,'first') - 1;
        sup2 = min(X2) + find(hau2sign==1,1,'last');
        
        
        %%%%% Measures %%%%%%
        msr_pi_max_x = sup1;
        msr_pi_max_y = smoothfit12(sup1);
        msr_hau_max_x = inf2;
        msr_hau_max_y = smoothfit22(inf2);
        
        ante12 = interp1(smoothfit12(inf1:sup1),inf1:sup1,1:101);
        ante22 = interp1(smoothfit22(inf2:sup2),inf2:sup2,1:101);
        area = ante22-ante12;
        
        msr_area_10 = sum(area(10:20));
        msr_area_20 = sum(area(20:30));
        msr_area_30 = sum(area(30:40));
        msr_area_40 = sum(area(40:50));
        msr_area_50 = sum(area(50:60));
        msr_area_60 = sum(area(60:70));
        msr_area_70 = sum(area(70:80));
        msr_area_80 = sum(area(80:90));
        msr_area_90 = sum(area(90:100));
        
        msr_energy = sum(Sreca(:));
        
        piarea = msr_pi_max_x-ante12;
        msr_picumarea = sum(piarea(20:find(ante12>0,1,'last')));
        
        PLOT = 1;
        if PLOT
            % Focus on the current signal
            figure(k)
            xlim([co2timeint(tinf)-1 co2timeint(tsup)+1])

            % Plot
            figure(2)
            subplot(121), plotmat(log(Sreca)); title('Raw') %Tint(trange),Fint(frange),
            hold on
            plot(inf1:sup1, smoothfit12(inf1:sup1).*(smoothfit12(inf1:sup1)>0), 'k', 'LineWidth',2)
            plot(inf2:sup2, smoothfit22(inf2:sup2).*(smoothfit22(inf2:sup2)>0), 'k', 'LineWidth',2)
            hold off
            subplot(122), plotmat(mpi+mhau1+mhau2); title('The matrix to treat') % plotmat(new) % Tint(trange),Fint(frange),
            hold on
            plot(inf1:sup1, smoothfit12(inf1:sup1).*(smoothfit12(inf1:sup1)>0), 'r', 'LineWidth',2)
            plot(inf2:sup2, smoothfit22(inf2:sup2).*(smoothfit22(inf2:sup2)>0), 'r', 'LineWidth',2)
            plot(msr_pi_max_x,msr_pi_max_y,'*r')
            plot(msr_hau_max_x,msr_hau_max_y,'*r')
            hold off
            %subplot(133), plotseg(seg{i},0,1,0,0) % plotmat(newdn)

            %figure(3)
            %plot(areasum,'-r'), hold on, %xlim([30 52]), ylim([400 1000])

            % Press key to continue
            a = 1;
            while a
                a = ~waitforbuttonpress;
            end

            %figure(3)
            %plot(areasum,'b'), hold on, %xlim([30 52]), ylim([400 1000])
        end
    catch err
        disp(err)
        disp('Nous maîtrisons la situation !')
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


