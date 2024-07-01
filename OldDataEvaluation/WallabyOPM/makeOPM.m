% function OPM = makeOPM(maps)
% OPM = reshape(mean(maps.*exp(2i*pi/size(maps,1)*reshape([1:size(maps,1)],size(maps,1),1,1)),1),size(maps,2),size(maps,3));
% end

function OPM = makeOPM(maps)
OPM = reshape(mean(maps.*exp(2i*pi/size(maps,3)*reshape([1:size(maps,3)],1,1,size(maps,3))),3),size(maps,1),size(maps,2));
end
