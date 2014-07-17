function template = readtemplate(filename,size)
template = double(rgb2gray(imread(filename)));
template = flipdim(template,1);
template = imresize(template,size);
end