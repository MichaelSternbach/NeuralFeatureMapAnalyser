function P = calcProbSmaller(set,value)
P=sum(set<=value,'all')/length(set);
end