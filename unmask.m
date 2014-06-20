function mask = unmask(x, y, num, mask)
    mask(y, x:x + num) = 0;
end