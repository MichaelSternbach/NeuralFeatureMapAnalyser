function plotPwTracks3D(PinwheelTracks)
    for ii = 1:size(PinwheelTracks.label,1)
        x=PinwheelTracks.x(ii,:);
        y=PinwheelTracks.y(ii,:);
        z = 1:size(PinwheelTracks.label,2);

        plot3(x, y, z, 'LineWidth', 1);   % 'LineWidth' is optional for thicker lines
        hold on
    end

    grid on; 

end