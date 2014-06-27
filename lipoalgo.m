%%%%%%%%%%%%%%%%%%%%%%%
%% Initalization
%%%%%%%%%%%%%%%%%%%%%%%
tic, clear, clc, close all
lipoalgopaths;

global fs nfft ovlp T F

% Read .wav file
%[x, fs] = audioread('lipo.WAV');
[x, fs] = audioread('120119_071_mono1.wav',[round(2400*44100) round(2600*44100)]);

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
Safn = nlfilter(Saf,[20 2],f); %%%%%%%%%%% /!\ ARBITRARY CONST in pixel %%%%%%%%%%%

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
shftw = time2co(1.2); % Pi + Hau average lengths %%%%%% /!\ CONST %%%%%%

pas = 1; % en pixel
sumit = zeros(ceil((size(Safn,2)-shftw)/pas),1);
for i = pas:pas:(size(Safn,2)-shftw)
    sumit(i/pas) = sum(sum(abs(Safn(freq2co(LOWR):freq2co(HIGHR),i:i+shftw))));
end

% Peak detection
delta = 20;         %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%
[maxtab, ~] = peakdet(sumit, delta, 1:length(sumit));

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
% Treatment 
%%%%%%%%%%%%%%%%%%%%%%%

for i = 2:length(maxtab(:,1))
    disp(i)
    
    % Get position (time) of calls
    val = maxtab(i,1);
    [tinf, tsup] = findtco(val, shftw); %%%%%%%%%%% /!\ ARBITRARY CONST IN FUNCTION %%%%%%%%%%%
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
    
    figure(1)
    xlim([co2time(tinf)-3 co2time(tsup)+3])

    %%%%%%%%%%%%%%%%%%%%%%%
    % Measurments
    %%%%%%%%%%%%%%%%%%%%%%%

    % Plot
    figure(2)
    subplot(131), plotmat(T(trange),F(frange),log(Sreca)); title('Raw')
    subplot(132), plotmat(m); title('The matrix to treat')

    mask = ones(size(m));
    mnew = m;
    pocs = {};

    while ~isequal(sum(mnew(:)),0)
        
        [val, ind] = max(mnew(:));
        [I, J] = ind2sub(size(mnew),ind);
        %I = co2freq(I+freq2co(LOWR));
        %J = co2time(J)+co2time(tinf);

       % New piece of curve (poc)
        poc = [];
        
        % Tolérance au décalage horizontale / verticale
        maxh = 5;  %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%
        maxv = 15; %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%
        
        for vdir = -1:2:1 % direction vertical (+1 monte, -1 descend)
            % disp('Initialisation x, y')
            x = J;
            y = I;

            % VERTICAL
            ENDV = 0;
            vdec = 0;

            while ~ENDV && vdec < maxv
                % HORIZONTAL
                hdec = 0;
                ENDH = 0;

                if isequal(mnew(y,x),0)
                    [xout, yout, bool] = lookleft(x, y, mnew);
                    if bool
                        % disp('Funky left')
                        x = xout;
                        y = yout;
                        ENDH = 1;
                    else
                        %disp('Going right')
                        hdir = 1;
                    end
                else
                    %disp('Going left')
                    hdir = -1;
                end

                while ~ENDH && hdec < maxh
                    %disp('Searching for junction …')
                    % tant qu'on n'a pas rencontré une jonction (0, ~0), on se décale (5x max) 
                    hdec = hdec + 1; % on incrémente le décalage
                    [x, y, ENDH] = hcheckandgo(x, y, mnew, hdir); % on tente le décallage
                    % ENDH = 2 : on est au bord de l'image
                    % ENDH = 1 : on a trouvé une jonction
                    % ENDH = 0 : on continue à se décaller
                end

                if isequal(ENDH,1)
                    poc = [poc; x y]; % append (x, y) to the list 
                    [x2, ~] = checkco(x + 15, y, mnew);  %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%
                    mask(y, x:x2) = 0; % update mask
                    vdec = 0;
                elseif isequal(ENDH,0) % no junction
                    x = x - hdir*hdec;
                    mask(y,x) = 0;
                    vdec = vdec + 1;
                elseif isequal(ENDH,2) % x edges
                    mask(y,x) = 0;
                    vdec = vdec + 1;
                end
                
                mnew = mnew.*mask; % update mnew
                
                % vcheckandgo
                [~, ynew] = checkco(x, y + vdir, mnew); % se décalle verticalement
                if isequal(ynew, y)
                    ENDV = 1; % on a atteind un bord
                else
                    y = ynew;
                end
            end
        end
        % Enregistre la portion de courbe si + de 15 px
        if size(poc, 1) > 15  %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%
            pocs{end+1} = struct('data',poc);
        end
    end

    colors = {'.y','.g','.b','.m','.c','.k'};
%     out = zeros(size(m)*100);
    
    for j = 1:size(pocs,2)
        [xi,yi] = readpocs(pocs{j}.data);
        
        % Normal
%         yy = smooth(xi,yi,0.3,'rloess'); 
%         figure(2), subplot(133)
%         plot(xi,yi,colors{rem(j,6)+1}), xlim([0 78]),ylim([0 101]), hold on
%         sortpoc = sortrows(pocs{j}.data,1);
%         xyb = sortpoc(round(size(pocs{j}.data,1)/2),:);
%         plot(xyb(1),xyb(2),'*r')
%         plot(xi,yy,'r-','LineWidth',2)
        
        % Reverse
        [yy,ind] = sort(yi);
        xx = smooth(yi,xi,0.3,'rloess'); 
        figure(2), subplot(133)
        plot(xi,yi,colors{rem(j,6)+1}), xlim([0 78]), ylim([0 101]), hold on
        plot(xx(ind),yy,'r-','LineWidth',2)
        
%         [xx100, yy100] = checkco(round(xx(ind)*100), round(yy*100), out);
%         out(sub2ind(size(out),yy100,xx100)) = 1;
        fit = createFitLin(xi,yi);
        pocs{j}.crease = sign(fit.p1);
        %plot(fit)
    end
    hold off

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

    a = 1;
    while a
        a = ~waitforbuttonpress;
    end
end