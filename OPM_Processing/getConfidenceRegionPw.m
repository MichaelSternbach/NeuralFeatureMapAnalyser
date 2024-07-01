function SizesCI = getConfidenceRegionPw(pinwheel_stats,field_size_pix,Confidence)
    if nargin <3
        Confidence = 0.95;
    end
    SizesCI = zeros([0 getN_PW(pinwheel_stats)]);
    for i_pw = 1:getN_PW(pinwheel_stats)
        SizeCI = plotPinwheel(pinwheel_stats.x(i_pw,:),pinwheel_stats.y(i_pw,:),pinwheel_stats.probability(i_pw,:),field_size_pix,Confidence);
        SizesCI(i_pw)= SizeCI;
    end
end

% function plotPinwheelStats(pinwheel_stats,field_size_pix)
%     for i_pw = 1:getN_PW(pinwheel_stats)
%         plotPinwheel(pinwheel_stats.x(i_pw,:),pinwheel_stats.y(i_pw,:),pinwheel_stats.probability(i_pw,:),field_size_pix)
%     end
% end

function  SizeCI = plotPinwheel(PWx,PWy,ProbabilityPW,field_size_pix,Confidence)
    ProbabilityLimitPW = .0;
    if ProbabilityPW >= ProbabilityLimitPW
        %plotPosition(PWx(1),PWy(1),ProbabilityPW)
        [PWx,PWy] = getCIPwPos(PWx,PWy,Confidence);
        CI = plotConfidenceRegion(PWx,PWy,field_size_pix);
        SizeCI = sum(CI,'all');
    end
end

function [PWx,PWy] = getCIPwPos(PWx,PWy,Confidence)
    d2=(PWy-PWy(1)).^2+(PWx-PWx(1)).^2; d2sort=sort(rmmissing(d2));
    d2include=d2sort(1:round(length(d2sort)*Confidence));
    argin = ismember(d2,d2include);
    PWx = PWx(argin);
    PWy = PWy(argin);
end

function plotPosition(PWx,PWy,ProbabilityPW)
    plot(PWx,PWy,'.white')
    text(PWx,PWy,num2str(ProbabilityPW),'Color','white','FontSize',5)
end
function CI=plotConfidenceRegion(PWx,PWy,field_size_pix)
    CI = points_confidence_region(PWx,PWy,field_size_pix,'hull');
    contour(CI,[1 1],'white','linewidth',0.3)
    %plotCI(CI)
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