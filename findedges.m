function BW = findedges(Sreca, Srecaf)
% Compute edge detection with all methods available on log and not matrix
% and on fuzzy and not matrix. Then remove points of low probability
% according to the usual shape of the envelope and computes the mean of all
% matrix found.
%  usage : BW = findedges(mat, fuzzy_mat)
%  BW is the output matrix of the same size of the input matrix
% -------------------------------------------------------------------------
    % Sobel
    [BWv] = edge(Sreca,'sobel',0.5,'vertical');
    [BWh] = edge(Sreca,'sobel',0.5,'horizontal');
    [BWvl] = edge(log(Sreca),'sobel',0.5,'vertical');
    [BWhl] = edge(log(Sreca),'sobel',0.5,'horizontal');
    
    [BWvf] = edge(Srecaf,'sobel',0.5,'vertical');
    [BWhf] = edge(Srecaf,'sobel',0.5,'horizontal');
    [BWvlf] = edge(log(Srecaf),'sobel',0.5,'vertical');
    [BWhlf] = edge(log(Srecaf),'sobel',0.5,'horizontal');
    
    % Canny
    [BWc] = edge(Sreca,'canny');
    [BWcl] = edge(log(Sreca),'canny');
    
    [BWcf] = edge(Srecaf,'canny');
    [BWclf] = edge(log(Srecaf),'canny');
    
    % Prewitt
    [BWp] = edge(Sreca,'prewitt');
    [BWpl] = edge(log(Sreca),'prewitt');

    [BWpf] = edge(Srecaf,'prewitt');
    [BWplf] = edge(log(Srecaf),'prewitt');
    
    % Robert
    [BWr] = edge(Sreca,'roberts');
    [BWrl] = edge(log(Sreca),'roberts');

    [BWrf] = edge(Srecaf,'roberts');
    [BWrlf] = edge(log(Srecaf),'roberts');
    
    % Log
    [BWl] = edge(Sreca,'log');
    [BWll] = edge(log(Sreca),'log');

    [BWlf] = edge(Srecaf,'log');
    [BWllf] = edge(log(Srecaf),'log');

    
    m = ones(size(BWv));
    % usage => drawlineinmat(m, x1, y1, x2, y2)
    
    % Top right corner
    m = drawlineinmat(m, size(m,2)/2, size(m,1), size(m,2), size(m,1)/2);
    m = drawlineinmat(m, round(size(m,2)*3/4), size(m,1), size(m,2), round(size(m,1)*1/4));
    m = drawlineinmat(m, round(size(m,2)*2/3), size(m,1), size(m,2), round(size(m,1)*1/4));
    m = drawlineinmat(m, round(size(m,2)*2/3), size(m,1), size(m,2), round(size(m,1)*3/4));

    % Top left corner
    m = drawlineinmat(m, 1, round(size(m,1)*3/4),round(size(m,2)*1/3), size(m,1));
    m = drawlineinmat(m, 1, round(size(m,1)*1/4),round(size(m,2)*1/6), size(m,1));
    m = drawlineinmat(m, 1, round(size(m,1)*1/5),round(size(m,2)*1/6), size(m,1));
    
    % Vertical under curve
    m = drawlineinmat(m, round(size(m,2)*28/76), 1, round(size(m,2)*32/76), round(size(m,1)*72/110));
    m = drawlineinmat(m, round(size(m,2)*18/76), 1, round(size(m,2)*25/76), round(size(m,1)*63/110));
    
    m = drawlineinmat(m, round(size(m,2)*46/76), 1, round(size(m,2)*49/76), round(size(m,1)*53/110));
    m = drawlineinmat(m, round(size(m,2)*54/76), 1, round(size(m,2)*49/76), round(size(m,1)*53/110));
    
    BW = 0;
    BWs = {BWv, BWh, BWvl, BWhl,... 
        BWvf, BWhf, BWvlf, BWhlf,... 
        BWc, BWcl, BWcf, BWclf,... 
        BWp, BWpl, BWpf, BWplf,... 
        BWr, BWrl, BWrf, BWrlf,... 
        BWl, BWll, BWlf, BWllf};
        
    for i = 1:length(BWs)
        BWs{i} = BWs{i}.*m;
        BWs{i} = bwareaopen(BWs{i}, 35); %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%
        BW = BW + double(BWs{i});
    end

    BW = BW./(length(BWs)-4);
%     BW = (double(BWv) + double(BWh) + double(BWvl) + double(BWhl) ...
%         + double(BWvf) + double(BWhf) + double(BWvlf) + double(BWhlf)...
%         + double(BWc) + double(BWcl) + double(BWcf) + double(BWclf) ...
%         + double(BWp) + double(BWpl) + double(BWpf) + double(BWplf) ...
%         + double(BWr) + double(BWrl) + double(BWrf) + double(BWrlf) ...
%         + double(BWl) + double(BWll) + double(BWlf) + double(BWllf))./20;
      
      
      
      
      
      
      
      
      
end