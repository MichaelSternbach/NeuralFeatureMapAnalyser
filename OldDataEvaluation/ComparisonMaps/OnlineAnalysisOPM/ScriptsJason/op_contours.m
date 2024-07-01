function c = op_contours(op, orientations)
%OP_CONTOURS Calculate the contours of an OP map for pretty plotting.
%
% c = op_contours(op, orientations)
%
% Calculates the iso-orientation contours of the given OP map. The orientations
% to use should be given in degrees, not radians. Used for plotting only.
%
% Example:
% c = op_contours(op, [-90, -67.5, -45, -22.5, 0, 22.5, 45, 67.5]);
% figure;
% hold on;
% cmap = hsv(8);
% for k = 1:length(c)
%     plot(c{k}(1,:), c{k}(2,:), 'Color', cmap(k,:));
% end
% hold off;

op= exp(1i*op/256*2*pi);

orientations = 2*orientations*pi/180;
c = cell(length(orientations), 1);
for k = 1:length(orientations)
    thetas = wrapToPi(angle(op) - orientations(k));
    base_contours = contourc(thetas, [0 0]);
    correct_contours = [];
    for j = 1:size(base_contours, 2)
        x = base_contours(1,j);
        y = base_contours(2,j);
        if x == 0
            correct_contours = [correct_contours, [nan; nan]];
            continue;
        end
        t = thetas(round(y), round(x));
        if abs(abs(t) - pi) < pi/4
            correct_contours = [correct_contours, [nan; nan]];
            continue;
        end
        correct_contours = [correct_contours, [x; y]];
    end
    c{k} = correct_contours;
end