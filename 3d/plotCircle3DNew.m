function [points] = plotCircle3DNew(center,normal,radius,color, linewidthp, plotflag)
theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
if plotflag
    plot3(points(1,:),points(2,:),points(3,:),color,'linewidth',linewidthp);
end
end