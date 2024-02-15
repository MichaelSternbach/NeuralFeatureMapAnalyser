function plotAbsTile(Data, title,maxMap,minMap,ROI,PwInfo,z,color_contur,linewidth,Xlim,Ylim)
    if nargin <6
        PwInfo = 0;
    end
    if nargin <7
        z=0;
    end
    if nargin<8
        color_contur = 'w';
    end
    if nargin <9
        linewidth = 1;
    end
    if nargin<10
        Xlim=0;
    end
    if nargin<11
        Ylim=0;
    end
    ax = nexttile;
    plot_mapAbs(Data,title,maxMap,minMap,ROI,ax)
    if isstruct(PwInfo)
        hold on; plot(PwInfo.pinwheel_stats.x(:,1),PwInfo.pinwheel_stats.y(:,1),'o','Color',color_contur,'linewidth',linewidth)
    end
    if length(z)>1
        hold on
        contour(real(z),[0 0],'linewidth',linewidth,'Color',color_contur)
        hold on
        contour(imag(z),[0 0],'linewidth',linewidth,'Color',color_contur)
    end
    
    if length(Xlim)==2
        xlim(Xlim)
    end
    if length(Ylim)==2
        ylim(Ylim)
    end
end