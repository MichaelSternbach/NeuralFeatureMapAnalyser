function plotAbsTile(Data, title,maxMap,minMap,ROI,PwInfo)
    if nargin <6
        PwInfo = 0;
    end
    ax = nexttile;
    plot_mapAbs(Data,title,maxMap,minMap,ROI,ax)
    if isstruct(PwInfo)
        for ii = 1:size(PwInfo.PWxList,1)
            hold on
            if PwInfo.signlist(ii) > 0
                marker = '.white';
            else
                marker = '.black';
            end
            plot(PwInfo.PWxList(ii),PwInfo.PWyList(ii),marker)
        end
    end
end