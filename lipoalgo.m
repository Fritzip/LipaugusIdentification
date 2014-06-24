%%%%%%%%%%%%%%%%%%%%%%%
%% Initalization
%%%%%%%%%%%%%%%%%%%%%%%
tic, clear, clc, close all
lipoalgopaths;

global fs nfft ovlp T F

% Read .wav file
%[x, fs] = audioread('lipo.WAV');
[x, fs] = audioread('120119_071_mono3.WAV',[round(2400*44100) round(2600*44100)]);

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
    m = rmnoisepts(BWs,35); %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%
    
   
    figure(1)
    xlim([co2time(tinf)-3 co2time(tsup)+3])
    
    
    mod = zeros(size(m));
    mask = ones(size(m));


    %%%%%%%%%%%%%%%%%%%%%%%
    % Measurments
    %%%%%%%%%%%%%%%%%%%%%%%
    tic

    mask = ones(size(m));
    mnew = m;
    
    % Plot
    figure(2)
    subplot(131), plotmat(T(trange),F(frange),log(Sreca)); title('Raw')
    subplot(132), plotmat(m); title('The matrix to treat')

    pocs = {};

    while ~isequal(sum(mnew(:)),0)
        
        [val, ind] = max(mnew(:));
        [I, J] = ind2sub(size(mnew),ind);
        %I = co2freq(I+freq2co(LOWR));
        %J = co2time(J)+co2time(tinf);

        maxh = 5;
        maxv = 3;

        poc = []; % piece of curve
        for vdir = -1:2:1 % direction vertical (+1 monte, -1 descend)
            % Initialisation x, y
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
                        % Funky left
                        x = xout;
                        y = yout;
                        ENDH = 1;
                    else
                        % Going right
                        hdir = 1;
                    end
                else
                    % Going left
                    hdir = -1;
                end

                while ~ENDH && hdec < maxh
                    % Searching for junction …
                    % tant qu'on n'a pas rencontré une jonction (0, ~0), on se décale (5x max) 
                    hdec = hdec + 1; % on incrémente le décalage
                    [x, y, ENDH] = checkandgo(x, y, mnew, hdir); % on tente le décalalage
                    % ENDH = 1 : on a trouver une jonction
                    % ENDH = 0 : on continue à se décaller OU on est au bord de l'image
                end

                if isequal(ENDH,1)
                    % Found a junction
                    % on poursuit
                    poc = [poc; x y]; % append (x, y) to the list 
                    [x2, ~] = checkco(x + 15, y, mnew);
                    mask(y, x:x2) = 0; % update mask
                    %figure(2), subplot(133), plotmat(mnew), hold on, plot(x,y,'*g'), hold off 
                    vdec = 0;
                else
                    x = x - hdir*hdec;
                    mask(y,x) = 0;
                    vdec = vdec + 1;
                end
                mnew = mnew.*mask;
                [~, ynew] = checkco(x, y + vdir, mnew); % se décalle verticalement
                if isequal(ynew, y)
                    ENDV = 1; % on a atteind un bord
                else
                    y = ynew;
                end
            end
        end
        if size(poc, 1) > 15
            pocs{end+1} = poc;
        end
    end
    toc

    colors = {'*k','*g','*b','*r','*m','*y','*c'};
    for j = 1:size(pocs,2)
        figure(2), subplot(133), plot(pocs{j}(:,1),pocs{j}(:,2),colors{rem(j,7)+1}), ylim([0 101]), hold on
    end
    hold off


    a = 1;
    while a
        a = ~waitforbuttonpress;
    end
end