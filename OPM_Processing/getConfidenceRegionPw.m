function SizesCI = getConfidenceRegionPw(pinwheel_stats,field_size_pix,Confidence,plot_CI,Plot_CI_Size_Scale,plotProb,color,linewidth,method)
    if nargin <3
        Confidence = 0.95;
    end
    if nargin <4
        plot_CI = true;
    end
    if nargin <5
        Plot_CI_Size_Scale = false;
    end
    if ~isnumeric(Plot_CI_Size_Scale) && Plot_CI_Size_Scale
        Plot_CI_Size_Scale = 1;
    end
    if nargin <6
        plotProb= false;
    end
    if nargin <7
        color = 'white';
    end
    if nargin < 8
        linewidth = 0.3;
    end
    if nargin < 9
        method = 'gaussians';
    end
    plot(pinwheel_stats.x(:,1),pinwheel_stats.y(:,1),'.',Color='white',MarkerSize=0.5)
    offset = 1;
    SizesCI = zeros([0 getN_PW(pinwheel_stats)]);
    for i_pw = 1:getN_PW(pinwheel_stats)
        SizeCI = plotPinwheel(pinwheel_stats.x(i_pw,:),pinwheel_stats.y(i_pw,:),pinwheel_stats.probability(i_pw,:),field_size_pix,Confidence,plot_CI,Plot_CI_Size_Scale,color,linewidth,method);
        SizesCI(i_pw)= SizeCI;
        if plotProb
            text(pinwheel_stats.x(i_pw,1)+offset,pinwheel_stats.y(i_pw,1)+offset,num2str(round(pinwheel_stats.probability(i_pw),3)),Color=color)
        end
    end
end

% function plotPinwheelStats(pinwheel_stats,field_size_pix)
%     for i_pw = 1:getN_PW(pinwheel_stats)
%         plotPinwheel(pinwheel_stats.x(i_pw,:),pinwheel_stats.y(i_pw,:),pinwheel_stats.probability(i_pw,:),field_size_pix)
%     end
% end

function  SizeCI = plotPinwheel(PWx,PWy,ProbabilityPW,field_size_pix,Confidence,plot_CI,Plot_CI_Size_Scale,color,linewidth,method)
    if nargin <8
        color = 'white';
    end
    if nargin <9
        linewidth = 0.3;
    end
    if nargin <10
        method = 'hull';
    end
    ProbabilityLimitPW = .0;
    if ProbabilityPW >= ProbabilityLimitPW
        %plotPosition(PWx(1),PWy(1),ProbabilityPW)
        CI = plotConfidenceRegion(PWx,PWy,field_size_pix,plot_CI,color,linewidth,method,Confidence);
        SizeCI = sum(CI,'all');
        if isnumeric(Plot_CI_Size_Scale)
            text(PWx(1)+1,PWy(1)+1,num2str(round(sqrt(SizeCI)*Plot_CI_Size_Scale)),'Color',color,'FontSize',8)
        end
    end
end


function plotPosition(PWx,PWy,ProbabilityPW)
    plot(PWx,PWy,'.white')
    text(PWx,PWy,num2str(ProbabilityPW),'Color','white','FontSize',5)
end
function CI=plotConfidenceRegion(PWx,PWy,field_size_pix,plot_CI,color,linewidth,method,Confidence)
    if nargin <5
        color = 'white';
    end
    if nargin <6
        linewidth = 0.3;
    end
    CI = points_confidence_region(PWx,PWy,field_size_pix,method,Confidence);
    if plot_CI
        contour(CI,[1 1],color,'linewidth',linewidth)
        %plotCI(CI)
    end
        
end

function N_PW = getN_PW(pinwheel_stats)
    N_PW = size(pinwheel_stats.x,1);
end

function plotCI(CI)
    for i_x = 1:size(CI,1)
        for i_y = 1:size(CI,2)
            if CI(i_x,i_y)==1
                plot(i_x,i_y,'.white')
            end
        end
    end
end