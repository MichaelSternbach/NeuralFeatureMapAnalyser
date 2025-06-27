

function plotMapData(data_obj,BloodVesselImg,average_spacing_mm,filter,folder,mm,width_scale_pix)
    
    if nargin <4
        filter = false;
    end
    if nargin <5
        folder = '';
    end
    if nargin <6
        mm = .5;
    end
    if nargin <7
        width_scale_pix = 5;
    end
    
    
    if filter == true
        Ncolumns = 2;
    else
        Ncolumns = 1;
    end

    f = figure();
    t = tiledlayout(Ncolumns,3);
    f.Position = [100 100 1100 400];

    nexttile;
    PlotBloodVessels(BloodVesselImg,ones(size(BloodVesselImg)),1)
    hold on
    contour(data_obj.ROI,[1 1],'red')

    spacing_pix = mm * data_obj.info.pix_per_mm;
    hold on
    plot([width_scale_pix,width_scale_pix+spacing_pix],[width_scale_pix, width_scale_pix],'-red')
    hold on
    text(width_scale_pix+spacing_pix/2-8,width_scale_pix+4,[num2str(mm) ' mm'],'Color','red')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    title('Bloodvessel Map')

    ax2 = nexttile;
    z_input = data_obj.read_map();
    scale = (data_obj.info.pix_per_mm*average_spacing_mm).^-1;
    plot_map(z_input,data_obj.ROI,0,1)
    spacing_pix = 1/scale;
    hold on
    plot([width_scale_pix,width_scale_pix+spacing_pix],[width_scale_pix,width_scale_pix],'-white')
    hold on
    text(width_scale_pix+spacing_pix/2,width_scale_pix+4,'Λ','Color','white')
    title(ax2,'preferred orientation not filtered')

    ax3 = nexttile;
    Abs = abs(z_input);
    plot_mapAbs(Abs,'selectivity not filtered',max(Abs(data_obj.ROI),[],'all'),min(Abs(data_obj.ROI),[],'all'),data_obj.ROI,ax3)
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    
    if filter == true
        
        ax4 = nexttile;
        z_input = data_obj.filter_map(data_obj.read_map());
        scale = (data_obj.info.pix_per_mm*average_spacing_mm).^-1;
        plot_map(z_input,data_obj.ROI,0,1)
        spacing_pix = 1/scale;
        hold on
        plot([width_scale_pix,width_scale_pix+spacing_pix],[width_scale_pix,width_scale_pix],'-white')
        hold on
        text(width_scale_pix+spacing_pix/2,width_scale_pix+4,'Λ','Color','white')
        title(ax4,'preferred orientation filtered')

        ax5 = nexttile;
        Abs = abs(z_input);
        plot_mapAbs(Abs,'selectivity filtered',max(Abs(data_obj.ROI),[],'all'),min(Abs(data_obj.ROI),[],'all'),data_obj.ROI,ax5)
        set(gca,'xtick',[])
        set(gca,'ytick',[])
    end
    
    if filter
            CIFile = [folder 'MapData' data_obj.info.ID];
        else
            CIFile = [folder 'MapDataNotFiltered' data_obj.info.ID];
    end
    
    print(f,'-depsc', [CIFile '.eps'])
    savefig(f,[CIFile '.fig'])
end