function CI = getConfidenceMeanBS(data,nBoot,alpha)

    %% get Bootstrap samples and jackknife samples
    [bootstrapSamples, jackknifeSamples] = resample_bootstrap_jackknife(data, nBoot)
    
    %% calculate mean
    BS_mean = mean(bootstrapSamples,1);
    jack_mean = mean(jackknifeSamples,1);
    
    CI = bootstrap_ci(BS_mean,mean(data(:)),jack_mean,alpha);


end



