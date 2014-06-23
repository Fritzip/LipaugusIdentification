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

for i = 2:length(maxtab(:,1))
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
    m = rmnoisepts(BWs,35); %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%
    
    % Transform m to 3 vectors x,y (positions) and w the weighted vector
    %[x, y, w] = mat3vec(m);

    figure(1)
    xlim([co2time(tinf)-3 co2time(tsup)+3])
    
    % Plot
    figure(2)
    subplot(131), plotmat(T(trange),F(frange),log(Sreca)); title('Raw')
    subplot(132), plotmat(m); title('The matrix to treat'), hold on % T(trange),F(frange),
    
    
    
    
%     [x, y, w, h] = checkcorec(J-6,I-10,12,20,m);
%     rectangle('Position',[x, y, w, h],...
%             'LineWidth',2,...
%             'LineStyle','--',...
%             'EdgeColor','r')
        
    %m(m<0.15*max(m(:))) = 0;    
    
    mod = zeros(size(m));
    mask = ones(size(m));
    
    %subplot(133), plotmat(m); title('15% less'), hold on
    
    
    
%     srec = m(y:y+h,x:x+w);
%     [x, y, w] = mat3vec(srec);
%     
%     figure(3), plotmat(srec)
    %edg = findedges(srec);
    %BW = mergeBWs(edg,4);
    %disp(BW)
    %figure(4), plotmat(BW)

%     mr = m(:,:).^3;
%     %mr(mr<0.4*max(mr(:))) = 0;
%     [x, y, w] = mat3vec(mr);
%     figure(3), plotmat(mr), hold on,
%     [fitresult, gof] = createFitSmooth(x,y,w);
%     plot(fitresult), hold off
    





    %%%%%%%%%%%%%%%%%%%%%%%
    % Measurments
    %%%%%%%%%%%%%%%%%%%%%%%
    tic

    mask = ones(size(m));

    % Plot
    figure(2)
    subplot(131), plotmat(T(trange),F(frange),log(Sreca)); title('Raw')
    subplot(132), plotmat(m); title('The matrix to treat')
    %subplot(133), plotmat(m.*mask); title('With mask')

    pocs = {};

    while ~isequal(sum(m(:).*mask(:)),0)
        mnew = m.*mask;
        [val, ind] = max(mnew(:));
        [I, J] = ind2sub(size(mnew),ind);
        %I = co2freq(I+freq2co(LOWR));
        %J = co2time(J)+co2time(tinf);

        maxh = 4;
        maxv = 3;

        poc = []; % piece of curve
        for vdir = -1:2:1 % direction vertical (monte ou descend)
            % Initialisation x, y
            %disp('Init')
            x = J;
            y = I;

            % VERTICAL
            ENDV = 0;
            vdec = 0;

            while ~ENDV && vdec < maxv
                % HORIZONTAL
                hdec = 0;
                ENDH = 0;
                %sprintf('x = %d\ty = %d'
                if isequal(mnew(y,x),0)
                    [xout, yout, bool] = lookleft(x, y, mnew);
                    if bool
                        %disp('Funky left')
                        x = xout;
                        y = yout;
                        ENDH = 1;
                    else
                        %disp('Going right')
                        hdir = 1; % direction horizontale (gauche ou droite)
                    end
                else
                    %disp('Going left')
                    hdir = -1;
                end

                while ~ENDH && hdec < maxh
                    %disp('Searching for junction …')
                    % tant qu'on n'a pas rencontré une jonction (0, ~0), on se décale (5x max) 
                    hdec = hdec + 1; % on incrémente le décalage
                    [x, y, ENDH] = checkandgo(x, y, mnew, hdir); % on tente le décalalage
                    % ENDH = 1 : on a trouver une jonction
                    % ENDH = 0 : on continue à se décaller OU on est au bord de l'image
                end

                if isequal(ENDH,1)
                    %disp('Found a junction')
                    % on poursuit
                    poc = [poc; x y]; % append (x, y) to the list 
                    [x2, ~] = checkco(x + 10, y, mnew);
                    mask(y, x:x2) = 0; % update mask
                    %size(mask)
                    mnew = mnew.*mask;
                    %figure(2), subplot(133), plotmat(mnew), hold on, plot(x,y,'*g'), hold off 
                    %waitforbuttonpress;
                    vdec = 0;
                else
                    x = x - hdir*hdec;
                    vdec = vdec + 1;
                end
                [~, ynew] = checkco(x, y + vdir, mnew); % se décalle verticalement
                if isequal(ynew, y)
                    ENDV = 1; % on a atteind un bord
                else
                    y = ynew;
                end
            end
        end
        if size(poc, 1) > 10
            pocs{end+1} = poc;
        end
    end
    toc

    for i = 1:size(pocs,2)
        figure(2), subplot(133), plot(pocs{i}(:,1),pocs{i}(:,2),'*'), ylim([0 101]), hold on
    end
    hold off










    a = 1;
    while a
        a = ~waitforbuttonpress;
    end
end