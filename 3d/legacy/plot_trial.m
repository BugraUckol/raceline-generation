% t = 0:0.1:10;
% x = sin(t);
% s = sin(t).^2 + 2 * t;

close all
figure(8)
t = time_arr;
s = path.s_arr;
x = t;

% Bottom axis: x vs s
ax1 = axes;
plot(s, x, 'k', 'LineWidth', 2);
xlabel(ax1, 's [m]')
ylabel(ax1, 'x [m]')
title('x vs s with time as top axis')
grid on

% Top axis: t(s)
ax2 = axes;
plot(ax2, s, x, 'k', 'LineWidth', 2); 
ax2.Color = 'none';          
ax2.XAxisLocation = ['top'];     
ax2.YAxisLocation = 'right';  
ax2.YColor = 'none';   

% % Tick mapping
% xticks = get(ax1, 'XTick');
% t_interp = interp1(s, t, xticks);
% xticklabels = arrayfun(@(x) sprintf('%.1f', x), t_interp, 'UniformOutput', false);
% set(ax2, 'XTick', xticks, 'XTickLabel', xticklabels)

% Tick mapping
xtticks = interp1(t,s,t(1:100:end));
xticklabels = arrayfun(@(x) sprintf('%.1f', x), xtticks, 'UniformOutput', false);
set(ax2, 'XTick', xtticks, 'XTickLabel', t(1:100:end))

xlabel(ax2, 't(s) [s]')

% Ensure axes are aligned
linkaxes([ax1, ax2], 'xy')

% % Draw custom vertical grid lines on secondary axis
hold(ax2, 'on')
ylims = get(ax2, 'YLim');
for i = 1:length(xticks)
    xline(ax2, xtticks(i), '--r', 'LineWidth', 2);  % red dashed grid lines
end
hold(ax2, 'off')
