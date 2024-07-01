function maps = getMultiSourceMapESD(SourceList,simg,W,sigWin)
    N_Stimulus = 8;%size(W,1);
    maps = zeros([size(simg{1,1},1:2) N_Stimulus]);
    for Stimulus = 1: N_Stimulus
    	first = true;
        for Source = SourceList{Stimulus}
            disp([Stimulus Source])
            if first == true
            	maps(:,:,Stimulus) = (-simg{Stimulus,1}(:,:,Source))*mean(W{Stimulus,1}(sigWin,Source),1);
    	        first = false;
    	    else
    	    	maps(:,:,Stimulus) = maps(:,:,Stimulus)+(-simg{Stimulus,1}(:,:,Source))*mean(W{Stimulus,1}(sigWin,Source),1);
        end
    end
end
