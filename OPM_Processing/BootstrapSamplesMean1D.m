function samples_BS = BootstrapSamplesMean1D(samples,N_samplesBS,seed)
    n_samples = length(samples);
    samples_BS = zeros(1,N_samplesBS);
    rng(seed)
    for ii = 1:N_samplesBS
        samples_BS(ii) = mean(samples(randi(n_samples, n_samples, 1)),'all');
    end
end