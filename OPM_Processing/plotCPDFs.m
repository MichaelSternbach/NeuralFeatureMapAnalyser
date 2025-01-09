function plotCPDFs(data,label,LineMarker)
    if nargin<2
        label = '';
    end
    if nargin<3
        LineMarker = '-';
    end
    l = length(data);
    plot(sort(data),(1:l)/l,LineMarker,'DisplayName',label)
end