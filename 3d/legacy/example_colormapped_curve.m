figure(9)
clf
x = x_arr;
y = y_arr;
z = z_arr;
v = sqrt(u_arr.^2 + v_arr.^2 + w_arr.^2);
patch([x nan],[y nan],[z nan],[v nan], 'edgecolor', 'interp','LineWidth',5);
hcb = colorbar;colormap(jet);grid minor;
colorTitleHandle = get(hcb,'Title');
titleString = 'Total Velocity';
set(colorTitleHandle ,{'String','Rotation','Position'},{ titleString,90,[60 200]});
xlabel('x1'), ylabel('y1'), zlabel('Step Size [s]'), title('Step Size and Trajectory')
daspect([1,1,1])