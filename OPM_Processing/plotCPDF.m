function plotCPDF(data,label,ax)
    if nargin<2
        label = '';
    end
    if nargin<3
        ax = axes;
    end
    l = length(data);
    plot(ax,sort(data),(1:l)/l*100,'DisplayName',label)
end