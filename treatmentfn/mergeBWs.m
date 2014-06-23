function BW = mergeBWs(varargin)
%
%--------------------------------------------------------------------------
    BW = 0;
    if isequal(size(varargin,2),3)
        BWs = varargin{1};
        minpixgrp = varargin{2};
        mask = varargin{3};
    elseif isequal(size(varargin,2),2)
        BWs = varargin{1};
        minpixgrp = varargin{2};
        mask = ones(size(BWs{1}));
    end
    
    for i = 1:length(BWs)
        BWs{i} = BWs{i}.*mask;
        BWs{i} = bwareaopen(BWs{i}, minpixgrp); %%%%%%%%%%% /!\ ARBITRARY CONST %%%%%%%%%%%
        BW = BW + double(BWs{i});
    end

    BW = BW./(length(BWs)-2);
end