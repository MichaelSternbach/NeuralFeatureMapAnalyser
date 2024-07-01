function [ci,z0,acc] = bootstrap_ci(bootstat,stat,jackstat,alpha)
%{
Corrected and accelerated percentile bootstrap CI.
If bootstat is complex, circular statistics are applied to angles.

Equations from:
 Efron, B., & Hastie, T. (2016). Computer Age Statistical Inference:
        Algorithms, Evidence, and Data Science

Input:
 bootstat : [M,samples] array of bootstrap samples of statistic
 stat : [M,1] vector of measured statistic using all data. Each m of M is
       independent
 jackstat : [M,N] array of jackknife samples of statistic. N is the number
            of data points for each m.

Output:
 ci  : [M,lower,upper] array with confidence intervals for each M
 z0  : bias correction
 acc : acceleration constant

 Based on Matlab function bootci.m. Requires Statistics Toolbox.
 27 June 2016: chepe@nld.ds.mpg.de
%}

%% check if circular stats needed

if ~all(isreal(bootstat(:)))
    circ_flag = true;
    
    % convert to --> exp(1i*angle)
    bootstat = bootstat./abs(bootstat);
    
    % if no center given, calculate in one loop over columns the 
    % median angle for each row in bootstat
    if nargin<2
        
        % sum the differences from the angles of the first sample
        cir_center = bootstat(:,1);
        summed_difference = sum(abs(angle(bsxfun(@rdivide,bootstat,cir_center))),2);
        
        % go over the rest of the samples
        for ii=2:size(bootstat,2)
            
            % summed the differences from the samples of the current sample
            test_summed_difference = sum(abs(angle(bsxfun(@rdivide,bootstat,bootstat(:,ii)))),2);
            
            % check pixel indices that are closer to the median and exchange
            idx = find(test_summed_difference<summed_difference);
            summed_difference(idx)=test_summed_difference(idx);
            cir_center(idx) = bootstat(idx,ii);
        end
        
    else        
        cir_center = stat;
    end
    
    % convert to angles and center around stat
    bootstat = angle(bsxfun(@rdivide,bootstat,cir_center));
    stat = zeros(size(bootstat,1),1);
    
    % center jackknive samples around its mean
    if nargin >= 3
        % convert to --> exp(1i*angle)
        jackstat = jackstat./abs(jackstat);
        jackstat = angle(bsxfun(@rdivide,jackstat,sum(jackstat,2)));
    end
    
else
    circ_flag = false;
    
    % center jackknive samples around its mean
    if nargin >= 3
        jackstat = bsxfun(@minus,jackstat,mean(jackstat,2));
    end
end

%% calculate corrections and percentiles

% confidence level
if nargin < 4
    alpha = 0.05;
end

% Acceleration constant (Eq.11.40)
if nargin < 3
    acc = 0;
else
    acc = 1/6*sum(jackstat.^3,2)./(sum(jackstat.^2,2).^1.5);
    acc(all(jackstat==0,2)) = 0;
end

% bias-correction constant (Eq.11.31 - 11.32)
if nargin < 2
    z0 = 0;
else
    z0 = norminv(mean(bsxfun(@le,bootstat,stat),2));
end

% correct percentiles using constants (Eq.11.39)
z_alpha1 = norminv(alpha/2);
pct1 = 100*normcdf(z0 +(z0+z_alpha1)./(1-acc.*(z0+z_alpha1)));
pct1(z0==Inf) = 100;
pct1(z0==-Inf) = 0;

z_alpha2 = -z_alpha1;
pct2 = 100*normcdf(z0 +(z0+z_alpha2)./(1-acc.*(z0+z_alpha2)));
pct2(z0==Inf) = 100;
pct2(z0==-Inf) = 0;

%% get CI

M = size(bootstat,1);
ci = zeros(M,2);
if numel(pct1) == 1
    ci(:,1) = prctile(bootstat,pct2,2);
    ci(:,2) = prctile(bootstat,pct1,2);
else
    for m=1:M
        ci(m,1) = prctile(bootstat(m,:),pct2(m),2);
        ci(m,2) = prctile(bootstat(m,:),pct1(m),2);
    end
end
ci = sort(ci,2);

% if circular, shift back to initial position
if circ_flag
    ci = angle(bsxfun(@times,exp(1i*ci),cir_center));    
end

end
