% Read .wav file
clear, clc, close all
lipoalgopaths;
filename = '120119_071_mono3.wav';
info = audioinfo(filename);
step = 100;
cut = 2600;
measures = {};
%try
    %while cut*info.SampleRate < info.TotalSamples
for i = 1:5
    x = audioread(filename,[(cut-step)*info.SampleRate cut*info.SampleRate]);
    measures = [measures; lipoalgo(x,info.SampleRate, cut)];
    cut = cut + step;
end
%catch
    
%end
%%
disp('yipa')
csvwrite('measures.csv',measures)

D = pdist(measures,'cityblock');
L = linkage(D,'average');
dendrogram(L,0);

Y = measures(1:10,:);
D = pdist(Y,'cityblock');
L = linkage(D);
dendrogram(L,0);
c = cophenet(L,D)

