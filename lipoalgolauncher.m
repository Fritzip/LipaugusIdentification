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
    measures = [measures; lipoalgo(x,info.SampleRate, cut, i)];
    cut = cut + step;
end
%catch
    
%end
%%
disp('yipa')
csvwrite('measures.csv',measures(:,3:end))

D = pdist(cell2mat(measures(:,3:end)),'cityblock');
L = linkage(D,'average');
dendrogram(L,0,'labels',measures(:,2));

Y = measures(1:10,:);
D = pdist(Y,'cityblock');
L = linkage(D);
dendrogram(L,0);
c = cophenet(L,D)

{{'a' 'b' 'c' 'd' 'd' 'c' 'a' 'c' 'd' 'e' 'f' 'd' 'g' 'b' 'h' 'g' 'b'}... 
{'h' 'g' 'd' 'c' 'c' 'b' 'd' 'g' 'b' 'c' 'g' 'b' 'h' 'd' 'b' 'd' 'e'}...
{'c' 'Y' 'b' 'd' 'e' 'd' 'b' 'd' 'b' 'd' 'b' 'g' 'f' 'd' 'b' 'e' 'f' 'd' 'g' 'd' 'e' 'b' 'j' 'b' 'g'}...
{'a' 'c' 'j' 'b' 'a' 'g' 'j' 'b' 'e' 'd' 'j' 'b' 'g' 'a' 'j' 'b' 'e' 'b' 'j' 'a' 'b' 'j' 'g' 'a'}...
{'g' 'j' 'a' 'b' 'a' 'g' 'j' 'g' 'd' 'b' 'j' 'c' 'g' 'a' 'b' 'e' 'X' 'j' 'g' 'a' 'c' 'e' 'j' 'd' 'b'}}
















