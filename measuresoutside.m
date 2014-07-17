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

        au = (smoothfit1(1:0.1:30)).*(smoothfit1(1:0.1:30)>0);
        ie = (smoothfit2(1:0.1:40)).*(smoothfit2(1:0.1:40)>0);
        
        fitresults = createFitFourier2(areasum);
        measures = [measures; tinf double(uint8(truth{k}{ni}))...
                    fitresults.a0 fitresults.a1 fitresults.b1...
                    fitresults.a2 fitresults.b2 fitresults.w...
                    max(areasum) length(areasum) sum(Sreca(:))...
                    increasesize(areasum,size(m,1)-value) au' ie'];
                
        ni = ni+1;
    end