%% Launch lipoalgo
clear, clc, close all
lipoalgopaths;

filename = '120119_071_mono3.wav';
info = audioinfo(filename);
step = 100;
cut = 2600;
measures = [];
% while cut*info.SampleRate < info.TotalSamples
for i = 1:5
    x = audioread(filename,[(cut-step)*info.SampleRate cut*info.SampleRate]);
    measures = [measures; lipoalgo(x,info.SampleRate, cut, i)];
    cut = cut + step;
end
clearvars filename info step cut x i;
%% Write measures in .dat file
fileID = fopen('measures2.dat','w');
formatSpec = ['%d %s' repmat(' %f', [1,size(measures,2)-2]) '\n'];
[nrows,~] = size(measures);
for row = 1:nrows
    fprintf(fileID,formatSpec,measures(row,1),char(measures(row,2)), measures(row,3:end));
end
fclose(fileID);
clearvars formatSpec fileID nrows row;
%% Read measures in .dat file
% Open the text file.
fileID = fopen('measures.dat','r');
% Format string for each line of text:
formatSpec = '%f%s%f%f%f%f%f%f%f%[^\n\r]';

% Read columns of data according to format string.
data = textscan(fileID, formatSpec, 'Delimiter',' ');
classes = data(:,2);
classes = char(classes{1});
data = data(1:end-1);
data = dataset(data{1:end-1});
delcol = ~any(ismissing(data(:,3:end)),2);
delcol(classes=='X',:) = 0;
data = data(delcol,:);
classes = classes(delcol,:);
measures = double(data(:,3:end));
measures = cell2mat(cellfun(@(x)((x-mean(x))/std(x)), num2cell(measures,1),'UniformOutput',false));
% Close the text file.
fclose(fileID);
% Clear temporary variables
clearvars formatSpec fileID delcol;
 
%% Hierarchical Clustering
D = pdist(measures,'euclidean');
L = linkage(D,'ward');
dendrogram(L,0,'labels',classes); 

%% Kmeans
idx = kmeans(measures,5,'distance','city');
[silh3,h] = silhouette(measures,idx,'city');
set(get(gca,'Children'),'FaceColor',[.8 .8 1])
xlabel('Silhouette Value')
ylabel('Cluster')


% %% des tentatives de DFA
% clear, clc
%  load fisheriris
%  SL = meas(51:end,1);
%  SW = meas(51:end,2);
%  group = species(51:end);
%  h1 = gscatter(SL,SW,group, 'rb' , 'v^' ,[], 'off' );
%  set(h1, 'LineWidth' ,2)
%  legend( 'Fisher versicolor' , 'Fisher virginica' , ...
%         'Location' , 'NW' )
% 
% [X,Y] = meshgrid(linspace(4.5,8),linspace(2,4));
%  X = X(:);  Y = Y(:);
%  [C,err,P,logp,coeff] = classify([X Y],[SL SW], ...
%                                  group, 'Quadratic' );
%                              
%     %                       
% hold on ;
%  gscatter(X,Y,C, 'rb' , '.' ,1, 'off' );
%  K = coeff(1,2).const;
%  L = coeff(1,2).linear;
%  Q = coeff(1,2).quadratic;
%  
%  % Function to compute K + L*v + v'*Q*v for multiple vectors
%  % v=[x;y].  Accepts x and y as scalars or column vectors.
%  f = @(x,y) K + [x y]*L + sum(([x y]*Q) .* [x y], 2);
% 
%  h2 = ezplot(f,[4.5 8],[2 4]);
%  
%  set(h2, 'Color' , 'm' , 'LineWidth' ,2)
%  axis([4.5 8 2 4 ])
%  xlabel( 'Sepal Length' )
%  ylabel( 'Sepal Width' )
%  title( '{\bf Classification with Fisher Training Data}' )
%  
% %% R
% system('R.exe CMD BATCH test.R outfile');
% 
% %% DFA
% measures = cell2mat(cellfun(@(x)((x-mean(x))/std(x)), num2cell(measures,1),'UniformOutput',false));
% [C,err,P,logp,coeff] = classify(measures, measures, classes, 'Linear' );
% %%
% %[class,err,POSTERIOR,logp,coeff] = classify(measures,measures,classes);
% 
%  h1 = gscatter(measures(:,1), measures(:,3),classes, 'rb' , 'v^' ,[], 'off' );
%  set(h1, 'LineWidth' ,2)
% 
% 
% [X,Y] = meshgrid(linspace(-2,2),linspace(-2,2));
%  X = X(:);  Y = Y(:);
%  [C,err,P,logp,coeff] = classify([X Y],[measures(:,1) measures(:,3)], ...
%                                  classes, 'Quadratic' );
%                              
%     %                       
% hold on ;
%  gscatter(X,Y,C, 'rb' , '.' ,1, 'off' );
%  K = coeff(1,2).const;
%  L = coeff(1,2).linear;
%  Q = coeff(1,2).quadratic;
%  
%  % Function to compute K + L*v + v'*Q*v for multiple vectors
%  % v=[x;y].  Accepts x and y as scalars or column vectors.
%  f = @(x,y) K + [x y]*L + sum(([x y]*Q) .* [x y], 2);
% 
%  h2 = ezplot(f,[-2 2],[-2 2]);
%  
%  set(h2, 'Color' , 'm' , 'LineWidth' ,2)
%  axis([-2 2 -2 2])
%  xlabel( 'Sepal Length' )
%  ylabel( 'Sepal Width' )
%  title( '{\bf Classification with Fisher Training Data}' )






