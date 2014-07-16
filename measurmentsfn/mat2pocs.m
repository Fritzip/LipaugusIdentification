function [pocs, mout] = mat2pocs(m)
% Edges detector - get pieces of curve from matrix
    mask = ones(size(m));
    mout = zeros(size(m));
    mnew = m;
    pocs = {};
    
    while ~isequal(sum(mnew(:)),0)
        [~, ind] = max(mnew(:));
        [I, J] = ind2sub(size(mnew),ind);
        %I = co2freq(I+freq2co(LOWR));
        %J = co2time(J)+co2time(tinf);

       % New piece of curve (poc)
        poc = [];
        
        % Tol�rance au d�calage horizontale / verticale
        maxh = time2co(0.1);  %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%% 0.05 - 0.12
        maxv = freq2co(625); %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%
        
        for vdir = -1:2:1 % direction vertical (+1 monte, -1 descend)
            %disp('Initialisation x, y')
            x = J;
            y = I + double(vdir>0);

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
                        %disp('Funky left')
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
                    %disp('Searching for junction �')
                    % tant qu'on n'a pas rencontr� une jonction (0, ~0), on se d�cale (5x max) 
                    hdec = hdec + 1; % on incr�mente le d�calage
                    [x, y, ENDH] = hcheckandgo(x, y, mnew, hdir); % on tente le d�callage
                    % ENDH = 2 : on est au bord de l'image
                    % ENDH = 1 : on a trouv� une jonction
                    % ENDH = 0 : on continue � se d�caller
                end

                if isequal(ENDH,1)
                    poc = [poc; x y]; % append (x, y) to the list 
                    % Propagation average duration = 0.25 s
                    [xprop, ~] = checkco(x + time2co(0.25), y, mnew);  %%%%%%%%%%% /!\ CONST  %%%%%%%%%%%
                    mask(y, x:xprop) = 0; % update mask, remove the propagation behind the point found
                    mout(y, x) = m(y, x);
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
                [~, ynew] = checkco(x, y + vdir, mnew); % se d�calle verticalement
                if isequal(ynew, y)
                    ENDV = 1; % on a atteind un bord
                else
                    y = ynew;
                end
            end
        end
        % Enregistre la portion de courbe si + de 15 px
        if size(poc, 1) > freq2co(625)  %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%
            pocs{end+1} = struct('data', poc);
        end
    end
end