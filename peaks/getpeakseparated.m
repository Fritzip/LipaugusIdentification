function [pks, locs, over] = getpeakseparated(w,halfl,n,varargin)

% Works fine but very sloooow
%     forbideni = [];
%     pks = [];
%     locs = [];
%     for k = 1:n
%         m = 0;
%         while m < 1000
%             [C I] = max(w);
%             w(I) = 0;
%             if ~ismember(I,forbideni)
%                 forbideni = [forbideni I-halfl:I+halfl];
%                 pks = [pks C];
%                 locs = [locs I];
%                 break
%             end
%             m = m+1;
%             disp(m)
%         end
%     end


% Faster
    sort = 'descend';
    thresh = 0;
    
    if isequal(size(varargin,2),1) && ischar(varargin{1})
        sort = varargin{1};
    elseif isequal(size(varargin,2),1) && isnumeric(varargin{1})
        thresh = varargin{1};
    elseif isequal(size(varargin,2),2)
        if ischar(varargin{1}) && isnumeric(varargin{2})
            sort = varargin{1};
            thresh = varargin{2};
        elseif ischar(varargin{2}) && isnumeric(varargin{1})
            sort = varargin{2};
            thresh = varargin{1};
        end
    end

    
    [pks,locs] = findpeaks(w,'MINPEAKDISTANCE',2*halfl,'SORTSTR',sort);
    
    over = ones(length(pks),1);
    if thresh
        over(pks<thresh) = 0;
    end
    
    if length(locs) > n && ~isequal(n,0)
        pks = pks(1:n);
        locs = locs(1:n);
        over = over(1:n);
    end 
end