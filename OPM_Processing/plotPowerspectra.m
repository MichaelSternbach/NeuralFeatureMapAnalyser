function plotPowerspectra(power_profiles,color_rgb,name)

    for ii = 1:length(power_profiles)
        if ii==1
            plot(power_profiles{ii}.k_mm_inv,power_profiles{ii}.values_kspace, 'Color', [color_rgb 1])%,'DisplayName',name
        else
            plot(power_profiles{ii}.k_mm_inv,power_profiles{ii}.values_kspace, 'Color', [color_rgb 0.3])
        end
        hold on
    end
end