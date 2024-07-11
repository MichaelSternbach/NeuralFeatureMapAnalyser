function p = RankSumWidth(x,y)
    width_x = getWidth(x);
    width_y = getWidth(y);
    p = ranksum(width_x,width_y);
end

function width_x = getWidth(x)
    width_x = abs(x-mean(x));
end