set(gcf, 'color', 'none');
set(gca, 'color', 'none');
% exportgraphics(gcf,'transparent.svg',...   % since R2020a
% 'ContentType','vector',...
% 'BackgroundColor','none')
axis off;                    % remove axes, ticks, labels
set(gca, 'Visible', 'off');  % hide all decorations
title off
title('')
legend off