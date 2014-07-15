% Read .wav file
clear, clc, close all
filename = '120319_095_mono2.wav';
info = audioinfo(filename);
step = 100;
cut = 101;
m = 0;
%load('mastertemplate.mat')
while cut*info.SampleRate < info.TotalSamples
    x = audioread(filename,[(cut-step)*info.SampleRate cut*info.SampleRate]);
    m = m + createteamplate(x,m);
    cut = cut + step;
end
%save('mastertemplate.mat', 'm')