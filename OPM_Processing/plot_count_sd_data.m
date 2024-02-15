function plot_count_sd_data(areas, counts, plot_name)
% areas:     area in units of Lambda where the pinwheels are counted
% counts:    number of pinwheels in that area
% plot_name: name for the exported eps file 

% estimate count variance
no_bins = 100;
area_bins = logspace(log10(min(areas)),log10(max(areas)),no_bins+1); %% Logarithmically spaced area bins

sd = zeros(1,no_bins);
for jj = 1:no_bins
    rhos = counts((areas >= area_bins(jj)) & (areas < area_bins(jj+1))) ./ areas((areas >= area_bins(jj)) & (areas < area_bins(jj+1)));
    sd(jj) = std(rhos);
end

loglog(area_bins(1:end-1),sd,'markersize',12,'linewidth',2,'MarkerfaceColor','w','DisplayName',plot_name)

% xlim([0.3 20])
% ylim([0.1 5])
% ylabel('SD of pinwheel counts')
% xlabel('Region size [\Lambda]')
% set(gca,'xtick',[0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20])
% set(gca,'xticklabel',{'','','','','','','','1.0','','','','','','','','','10',''})
% set(gca,'ytick',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5])
% set(gca,'yticklabel',{'0.1','','','','','','','','','1.0','','','',''})
% set(gca,'Fontsize',25)
% text(8,3.5,'$$c(\frac{\rho}{A})^{\gamma}$$','Interpreter','latex','Fontsize',25)
% 
% % export to vector format file
% if nargin < 3
%     print('-depsc2', 'design_SD.eps');
% else
%     print('-depsc2', [plot_name,'.eps']);
% end

end