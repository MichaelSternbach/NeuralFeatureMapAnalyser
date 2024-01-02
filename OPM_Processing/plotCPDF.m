function plotCPDF(data,label,LineMarker,ax)
    if nargin<2
        label = '';
    end
    if nargin<3
        LineMarker = '-';
    end
    if nargin < 4
        ax = axes;
    end
    l = length(data);
    plot(ax,sort(data),(1:l)/l*100,LineMarker,'DisplayName',label)
end