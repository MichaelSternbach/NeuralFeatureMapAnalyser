function plot_nn_data_dist(d,edges, plot_name)
    % type:      d  , d++  , d+-
    % d:         measured distances in units of Lambda
    % plot_name: name for the exported eps file 
    if nargin <2
        edges = 0:0.05:1;
    end
    if nargin < 3
        plot_name = 'data'
    end

    % bin the measured distances
    [values,bins] = histcounts(d,edges);
    tmp = [bins(1:end-1);bins(2:end)];
    bins = tmp(:);
    tmp = [values;values];
    values = tmp(:);
    values = 100*values/sum(values)/2;

    % do plots

    plot(bins,values/max(values),'Linewidth',3,'DisplayName',plot_name)

    % xlim([0 1])
    % ylim([0 4.5])
    % xlabel('NN distance [\Lambda]')
    % ylabel('Freq. (normalized)')
    % set(gca,'xtick',[0 0.5 1])
    % set(gca,'ytick',[0 1 2 3 4])
    % set(gca,'yticklabel',{'0','','2','','4'})
    % text(0.05,4.1,type,'FontSize',40)
    % box off
    % set(gca,'Fontsize',25)
    % set(gca,'linewidth',1.5)
    % 
    % % export to vector format file
    % % if nargin < 3
    % %     print('-depsc2', 'design_nndist.eps');
    % % else
    % %     print('-depsc2', [plot_name,'.eps']);
    % % end

end