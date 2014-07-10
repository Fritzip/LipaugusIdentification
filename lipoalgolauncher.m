% Read .wav file
clear, clc, close all
lipoalgopaths;
filename = '120119_071_mono3.wav';
info = audioinfo(filename);
step = 100;
cut = 2600;
measures = [];
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
fileID = fopen('measures.dat','w');
formatSpec = ['%d %s' repmat(' %f', [1,size(measures,2)-2]) '\n'];
[nrows,ncols] = size(measures);
for row = 1:nrows
    fprintf(fileID,formatSpec,measures(row,1),char(measures(row,2)), measures(row,3:end));
end
fclose(fileID);

%csvwrite('measures.csv',measures)
%%
disp('yipa')
csvwrite('measures.csv',measures(:,3:end))
figure(6)
D = pdist(cell2mat(measures(:,3:end)));
L = linkage(D,'average');
dendrogram(L,0,'labels',measures(:,2));
ylim([0 3000])

c = cophenet(L,D)

{{'a' 'b' 'c' 'd' 'd' 'c' 'a' 'c' 'd' 'e' 'f' 'd' 'g' 'b' 'h' 'g' 'b'}... 
{'h' 'g' 'd' 'c' 'c' 'b' 'd' 'g' 'b' 'c' 'g' 'b' 'h' 'd' 'b' 'd' 'e'}...
{'c' 'Y' 'b' 'd' 'e' 'd' 'b' 'd' 'b' 'd' 'b' 'g' 'f' 'd' 'b' 'e' 'f' 'd' 'g' 'd' 'e' 'b' 'j' 'b' 'g'}...
{'a' 'c' 'j' 'b' 'a' 'g' 'j' 'b' 'e' 'd' 'j' 'b' 'g' 'a' 'j' 'b' 'e' 'b' 'j' 'a' 'b' 'j' 'g' 'a'}...
{'g' 'j' 'a' 'b' 'a' 'g' 'j' 'g' 'd' 'b' 'j' 'c' 'g' 'a' 'b' 'e' 'X' 'j' 'g' 'a' 'c' 'e' 'j' 'd' 'b'}}

double(uint8(truth{k}{ni}))

    [  56]    'a'    [   594.4568]    [ -748.1490]    [   614.5310]    [   132.4969]    [   45.1760]    [     0.0343]    [1.5835e+03]    [78]    [2.7507e+03]    [1x86 double]
    [ 172]    'b'    [   629.0332]    [ -735.0811]    [   388.0989]    [    91.7981]    [   82.6766]    [     0.0365]    [      1475]    [82]    [2.0075e+04]    [1x86 double]
    [ 510]    'c'    [   672.9387]    [ -790.5286]    [   370.5244]    [    96.0705]    [   99.3865]    [     0.0443]    [      1553]    [64]    [1.1228e+03]    [1x86 double]
    [ 647]    'd'    [   671.0891]    [ -801.4832]    [   455.4178]    [   113.5308]    [   91.2209]    [     0.0339]    [      1602]    [83]    [1.8374e+03]    [1x86 double]
    [1839]    'd'    [   566.4024]    [ -694.2998]    [   472.0306]    [   110.4775]    [   68.3372]    [     0.0328]    [      1408]    [83]    [2.3061e+03]    [1x86 double]
    [1935]    'c'    [   504.6463]    [ -640.0977]    [   504.7899]    [   113.2264]    [   43.5533]    [     0.0401]    [      1299]    [62]    [1.1207e+03]    [1x86 double]
    [2681]    'a'    [   549.6271]    [ -707.0091]    [   590.0693]    [   135.1535]    [   46.3757]    [     0.0347]    [      1492]    [78]    [2.2961e+03]    [1x86 double]
    [3225]    'c'    [   532.1432]    [ -656.5399]    [   438.2204]    [   104.0422]    [   51.4169]    [     0.0424]    [      1328]    [63]    [1.0393e+03]    [1x86 double]
    [3337]    'd'    [   153.3040]    [ -200.4614]    [   229.0484]    [    39.2401]    [   12.4467]    [     0.0351]    [       450]    [73]    [1.7604e+04]    [1x86 double]
    [3444]    'e'    [   342.4896]    [ -419.6533]    [   268.9178]    [    63.3255]    [   50.4104]    [     0.0378]    [       832]    [71]    [1.8411e+03]    [1x86 double]
    [3593]    'f'    [   502.9482]    [ -634.6656]    [   488.7068]    [   111.5500]    [   53.3084]    [     0.0365]    [      1319]    [77]    [2.0983e+03]    [1x86 double]
    [4127]    'd'    [   650.0922]    [ -768.5091]    [   419.7816]    [   102.4436]    [   89.8379]    [     0.0342]    [      1540]    [84]    [2.3068e+03]    [1x86 double]
    [4303]    'g'    [   604.9539]    [ -719.4434]    [   423.4578]    [    97.1454]    [   77.5117]    [     0.0381]    [      1452]    [76]    [1.0150e+03]    [1x86 double]
    [4406]    'b'    [   720.6070]    [ -849.8106]    [   458.7064]    [   112.2875]    [   97.3770]    [     0.0349]    [      1706]    [84]    [1.7542e+04]    [1x86 double]
    [5120]    'h'    [   565.2553]    [ -686.9178]    [   427.0409]    [   101.7796]    [   70.5296]    [     0.0386]    [      1394]    [74]    [1.2057e+03]    [1x86 double]
    [5236]    'g'    [   167.9859]    [ -244.0830]    [   373.7255]    [    65.3225]    [  -12.2832]    [     0.0353]    [       622]    [79]    [2.8806e+03]    [1x86 double]
    [5366]    'b'    [   663.0181]    [ -767.8899]    [   379.8795]    [    90.2821]    [   95.7296]    [     0.0354]    [      1528]    [83]    [1.7298e+04]    [1x86 double]
    
    

    [ 338]    'g'    [   690.5281]    [ -818.4801]    [   412.4307]    [   106.1400]    [  108.8155]    [     0.0388]    [      1628]    [76]    [2.8162e+03]    [1x86 double]
    [ 492]    'j'    [   526.2589]    [ -652.7708]    [   487.6344]    [   110.3647]    [   52.3529]    [     0.0336]    [      1356]    [85]    [9.4862e+03]    [1x86 double]
    [ 661]    'a'    [   466.2789]    [ -607.9318]    [   534.7037]    [   120.2993]    [   39.3170]    [     0.0352]    [      1291]    [78]    [2.4231e+03]    [1x86 double]
    [1291]    'b'    [   673.3506]    [ -780.7223]    [   383.4363]    [    92.3401]    [   98.4965]    [     0.0355]    [      1551]    [82]    [1.9725e+04]    [1x86 double]
    [1546]    'a'    [   514.2733]    [ -669.3354]    [   544.1057]    [   132.6989]    [   50.8208]    [     0.0357]    [      1400]    [78]    [2.4320e+03]    [1x86 double]
    [1710]    'g'    [   656.2977]    [ -707.1432]    [   215.6070]    [    40.7758]    [  107.3532]    [     0.0415]    [1.3995e+03]    [76]    [  818.7772]    [1x86 double]
    [1889]    'j'    [   604.1183]    [ -755.1791]    [   555.7979]    [   131.1299]    [   63.4646]    [     0.0330]    [      1565]    [85]    [1.0490e+04]    [1x86 double]
    [2100]    'g'    [   652.5111]    [ -801.5001]    [   514.6388]    [   124.8914]    [   81.1640]    [     0.0371]    [      1626]    [76]    [2.3932e+03]    [1x86 double]
    [2252]    'd'    [   679.9861]    [ -809.9406]    [   507.1379]    [   114.2103]    [   85.0075]    [     0.0323]    [      1629]    [84]    [2.3450e+03]    [1x86 double]
    [2560]    'b'    [   788.4562]    [ -919.8532]    [   475.9867]    [   114.4807]    [  109.7653]    [     0.0345]    [      1844]    [85]    [1.9320e+04]    [1x86 double]
    [2837]    'j'    [   606.8583]    [ -780.6938]    [   656.9143]    [   150.6219]    [   50.7816]    [     0.0322]    [      1653]    [85]    [1.0356e+04]    [1x86 double]
    [3103]    'c'    [   549.3823]    [ -669.5947]    [   411.4816]    [   100.3464]    [   65.4207]    [     0.0427]    [      1333]    [63]    [1.4193e+03]    [1x86 double]
    [3345]    'g'    [   574.9565]    [ -714.6444]    [   535.9552]    [   117.8850]    [   60.2535]    [     0.0361]    [      1478]    [76]    [2.0795e+03]    [1x86 double]
    [3587]    'a'    [   418.9079]    [ -569.4967]    [   599.9871]    [   128.2019]    [   18.2082]    [     0.0346]    [      1263]    [78]    [2.2665e+03]    [1x86 double]
    [3758]    'b'    [   743.6582]    [ -874.3395]    [   477.4925]    [   113.7221]    [  100.6597]    [     0.0346]    [      1759]    [85]    [1.8395e+04]    [1x86 double]
    [3926]    'e'    [   540.4850]    [ -692.0301]    [   578.2139]    [   130.9107]    [   43.7724]    [     0.0339]    [      1459]    [79]    [2.2733e+03]    [1x86 double]
    [3976]    'X'    [-4.3457e+11]    [5.7943e+11]    [-1.1059e+09]    [-1.4486e+11]    [5.5288e+08]    [-1.3595e-04]    [      1236]    [76]    [2.6634e+03]    [1x86 double]
    [4177]    'j'    [   674.3181]    [ -875.0037]    [   876.2500]    [   175.7147]    [   22.3479]    [     0.0299]    [      1939]    [86]    [1.0188e+04]    [1x86 double]
    [4587]    'g'    [   752.8792]    [ -934.7418]    [   692.8897]    [   155.7930]    [   77.1358]    [     0.0344]    [      1918]    [77]    [2.0386e+03]    [1x86 double]
    [4876]    'a'    [   466.0425]    [ -595.2024]    [   458.0454]    [   109.4247]    [   51.5571]    [     0.0372]    [      1233]    [77]    [2.3417e+03]    [1x86 double]
    [5073]    'c'    [   500.7434]    [ -633.2222]    [   500.3219]    [   110.2313]    [   41.8835]    [     0.0404]    [      1293]    [62]    [1.3635e+03]    [1x86 double]
    [5154]    'e'    [   617.9428]    [ -720.6584]    [   348.2716]    [    86.3869]    [   88.3809]    [     0.0373]    [      1433]    [79]    [1.9590e+03]    [1x86 double]
    [5305]    'j'    [   586.1510]    [ -750.6974]    [   694.2461]    [   142.9498]    [   33.6681]    [     0.0311]    [1.6275e+03]    [85]    [9.8360e+03]    [1x86 double]
    [5395]    'd'    [   692.3955]    [ -796.1925]    [   349.5642]    [    79.2813]    [   96.3546]    [     0.0350]    [1.5745e+03]    [84]    [3.0556e+03]    [1x86 double]
    [5515]    'b'














