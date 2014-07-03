function BW = rmnoisepts(BWs,minpixgrp)
% Remove points of low probability
% according to the usual shape of the envelope and computes the mean of all
% matrix found.
% -------------------------------------------------------------------------
    mask = ones(size(BWs{1}));
    % usage => drawlineinmat(m, x1, y1, x2, y2)
    
    % Top right corner
    mask = drawlineinmat(mask, size(mask,2)/2, size(mask,1), size(mask,2), size(mask,1)/2);
    mask = drawlineinmat(mask, round(size(mask,2)*3/4), size(mask,1), size(mask,2), round(size(mask,1)*1/4));
    mask = drawlineinmat(mask, round(size(mask,2)*2/3), size(mask,1), size(mask,2), round(size(mask,1)*1/4));
    mask = drawlineinmat(mask, round(size(mask,2)*2/3), size(mask,1), size(mask,2), round(size(mask,1)*3/4));

    % Top left corner
    mask = drawlineinmat(mask, 1, round(size(mask,1)*3/4),round(size(mask,2)*1/3), size(mask,1));
    mask = drawlineinmat(mask, 1, round(size(mask,1)*1/4),round(size(mask,2)*1/6), size(mask,1));
    mask = drawlineinmat(mask, 1, round(size(mask,1)*1/5),round(size(mask,2)*1/6), size(mask,1));
    
    % Vertical under curve
    mask = drawlineinmat(mask, round(size(mask,2)*28/76), 1, round(size(mask,2)*32/76), round(size(mask,1)*72/110));
    mask = drawlineinmat(mask, round(size(mask,2)*18/76), 1, round(size(mask,2)*25/76), round(size(mask,1)*63/110));
    
    mask = drawlineinmat(mask, round(size(mask,2)*46/76), 1, round(size(mask,2)*49/76), round(size(mask,1)*53/110));
    mask = drawlineinmat(mask, round(size(mask,2)*54/76), 1, round(size(mask,2)*49/76), round(size(mask,1)*53/110));
    
    % Horizontal criquet
%     mask = drawlineinmat(mask, 1, round(size(mask,1)*39/110), size(mask,2), round(size(mask,1)*39/110));
%     mask = drawlineinmat(mask, 1, round(size(mask,1)*47/110), size(mask,2), round(size(mask,1)*47/110));
%     mask = drawlineinmat(mask, 1, round(size(mask,1)*43/110), size(mask,2), round(size(mask,1)*43/110));
%     
    % 
    BW = mergeBWs(BWs,minpixgrp,mask);
end