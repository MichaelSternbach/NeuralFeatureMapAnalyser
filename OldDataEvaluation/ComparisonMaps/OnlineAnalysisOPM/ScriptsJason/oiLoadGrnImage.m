function grnImg = oiLoadGrnImage(filename),

  grnImg = single(imread(filename));

  mn = min(grnImg(:));
  mx = max(grnImg(:));

  grnImg = 255.*(grnImg-mn)./(mx-mn);