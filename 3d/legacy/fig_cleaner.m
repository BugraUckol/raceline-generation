set(gcf, 'color', 'none');
set(gca, 'color', 'none');
% exportgraphics(gcf,'empty_tube.eps',...   % since R2020a
% 'ContentType','vector',...
% 'BackgroundColor','none')
axis off;                    % remove axes, ticks, labels
set(gca, 'Visible', 'off');  % hide all decorations
title off
title('')
legend off
% exportgraphics(gcf,'empty_tube.eps',...   % since R2020a
% 'ContentType','vector',...
% 'BackgroundColor','none')

exportgraphics(gcf, 'tube_and_velocity_profile.png', 'Resolution',500)