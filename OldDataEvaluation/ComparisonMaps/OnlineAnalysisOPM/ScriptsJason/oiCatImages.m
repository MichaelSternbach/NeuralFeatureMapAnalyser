% concatanate image frames, one sequence/trial per "row"

[numConds, numTrials] = size(images);

[m,n,k] = size(images{1,1});

img = NaN([m*numTrials,n*k,numConds]);

for i = 1:numConds,
  for j = 1:numTrials,
    img((j-1)*m+1:j*m,:,i) = reshape(images{i,j},m,[]);
  end
end
