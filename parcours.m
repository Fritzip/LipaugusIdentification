plotmat(m), hold on
while ~isequal(sum(m.*mask),0)
    [val, ind] = max(m(:).*mask(:));
    [I, J] = ind2sub(size(m),ind);
    %I = co2freq(I+freq2co(LOWR));
    %J = co2time(J)+co2time(tinf);

    plot(J,I,'*g')
    hold on

    maxh = 5;
    maxv = 5;

    poc = []; % piece of curve
    for vdir = -1:2:1 % direction vertical (monte ou descend)
        % Initialisation x, y
        x = J;
        y = I;

        % VERTICAL
        ENDV = 0;
        vdec = 0;
        while ~ENDV && vdec <= maxv
            if isequal(m(y,x),0)
                hdir = 1; % direction horizontale (gauche ou droite)
            else
                hdir = -1;
            end

            % HORIZONTAL
            hdec = 0;
            ENDH = 0;
            while ~ENDH && hdec <= maxh
                % tant qu'on n'a pas rencontré une jonction (0, ~0), on se décale (5x max) 
                hdec = hdec + 1; % on incrémente le décalage
                [x, y, ENDH] = checkandgo(x, y, m, hdir); % on tente le décalalage
                % ENDH = 1 : on a trouver une jonction
                % ENDH = 0 : on continue à se décaller OU on est au bord de l'image
            end

            if isequal(ENDH,1)
                % on poursuit
                poc = [poc; x y]; % append (x, y) to the list 
                mask(y, x:x + 5) = 0; % update mask
                vdec = 0;
                [~, ynew] = checkco(x, y + vdir, m); % se décalle verticalement
                if isequal(ynew,y)
                    ENDV = 1; % on a atteind un bord
                else
                    y = ynew;
                end
            elseif isequal(ENDH,0)
                vdec = vdec + 1;
            end
        end
        % SAUVER POC > n pts
    end
end