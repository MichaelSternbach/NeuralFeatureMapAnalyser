function [pre,Order] = getOrder(Number)
 if Number ==0
     Order = 0;
 else
     Order = floor(log10(Number));
 end
 pre = Number/10^Order;
end