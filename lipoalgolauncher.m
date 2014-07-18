%% Launch lipoalgo
clear, clc, close all
lipoalgopaths;
filename = '120119_071_mono3.wav';
info = audioinfo(filename);
step = 100;
cut = 2600;
measures = [];
% while cut*info.SampleRate < info.TotalSamples
for i = 1:5 % »»»»» à augmenter «««««
    x = audioread(filename,[(cut-step)*info.SampleRate cut*info.SampleRate]);
    measures = [measures; lipoalgo(x,info.SampleRate, cut, i)];
    cut = cut + step;
end

%% Write measures in .dat file
fileID = fopen('.\output\measures3.dat','w');
formatSpec = ['%d %s' repmat(' %f', [1,size(measures,2)-2]) '\n'];
[nrows,ncols] = size(measures);
for row = 1:nrows
    fprintf(fileID,formatSpec,measures(row,1),char(measures(row,2)), measures(row,3:end));
end
fclose(fileID);

%% Hierarchical Clustering
%csvwrite('measures.csv',measures(:,3:end))
%measures(isnan(measures)) = 0;
D = pdist(measuresdfaonaeghf,'cityblock'); % measures(:,3:end)
L = linkage(D,'average');
dendrogram(L,0,'labels',char(dataV2')); % char(measures(:,2))
%c = cophenet(L,D)

%% Kmeans
idx = kmeans(mesuresaftlda2,9,'distance','city');
[silh3,h] = silhouette(mesuresaftlda2,idx,'city');
set(get(gca,'Children'),'FaceColor',[.8 .8 1])
xlabel('Silhouette Value')
ylabel('Cluster')







