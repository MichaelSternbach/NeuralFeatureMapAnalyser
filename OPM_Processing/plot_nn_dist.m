function plot_nn_dist(type, d, plot_name)
% type:      d  , d++  , d+-
% d:         measured distances in units of Lambda
% plot_name: name for the exported eps file 

% common design line
switch type
    case 'd'
        n = 1.2;
        xo = 0.48;
        b = 0.047;
    case 'd++'
        n = 4.5;
        xo = 0.5;
        b = 0.05;
    case 'd+-'
        n = 1.2;
        xo = 0.46;
        b = 0.08;
    otherwise 
        disp('Unknown type')
        return
end
func = @(x) x.^n./(1+exp((x-xo)/b));

% bin the measured distances
[values,bins] = histcounts(d,0:0.05:1);
tmp = [bins(1:end-1);bins(2:end)];
bins = tmp(:);
tmp = [values;values];
values = tmp(:);
values = 100*values/sum(values)/2;

% do plots
figure

distance = 0:0.01:1;

plot(distance,max(values)*func(distance)/max(func(distance)),'--k','Linewidth',2)
hold on
plot(bins,values,'k','Linewidth',3)

xlim([0 1])
ylim([0 4.5])
xlabel('NN distance [\Lambda]')
ylabel('Freq. (normalized)')
set(gca,'xtick',[0 0.5 1])
set(gca,'ytick',[0 1 2 3 4])
set(gca,'yticklabel',{'0','','2','','4'})
text(0.05,4.1,type,'FontSize',40)
box off
set(gca,'Fontsize',25)
set(gca,'linewidth',1.5)

% export to vector format file
% if nargin < 3
%     print('-depsc2', 'design_nndist.eps');
% else
%     print('-depsc2', [plot_name,'.eps']);
% end

end