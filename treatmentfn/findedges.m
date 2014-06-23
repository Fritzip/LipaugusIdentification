function BWs = findedges(mat)
%  Compute edge detection with all methods available on log and not matrix
%  usage : BWs = findedges(mat)
%  BWs is a cell array containing the outputs matrix of the same size of 
%  the input matrix
% -------------------------------------------------------------------------
    % Sobel
    [BWv] = edge(mat,'sobel',0.5,'vertical');
    [BWh] = edge(mat,'sobel',0.5,'horizontal');
    %[BWvl] = edge(log(mat),'sobel',0.5,'vertical');
    %[BWhl] = edge(log(mat),'sobel',0.5,'horizontal');
    
    % Canny
    [BWc] = edge(mat,'canny');
    %[BWcl] = edge(log(mat),'canny');
    
    % Prewitt
    [BWp] = edge(mat,'prewitt');
    %[BWpl] = edge(log(mat),'prewitt');

    % Robert
    [BWr] = edge(mat,'roberts');
    %[BWrl] = edge(log(mat),'roberts');
    
    % Log
    [BWl] = edge(mat,'log');
    %[BWll] = edge(log(mat),'log');

%     BWs = {BWv, BWh, BWvl, BWhl,... 
%         BWc, BWcl,... 
%         BWp, BWpl,... 
%         BWr, BWrl,... 
%         BWl, BWll};
    BWs = {BWv, BWh,... 
        BWc,... 
        BWp,... 
        BWr,... 
        BWl};
end