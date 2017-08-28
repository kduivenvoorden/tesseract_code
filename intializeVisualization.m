function intializeVisualization(L,f)
%undocumented file used for making pictures

global rs_global
rs = rs_global;

if length(rs)==4
    rs = rs(1:3);
end

boundaries = [(L-1)*ones(2,1);(f*(L-2)+1)]+rs;


f1 = figure('WindowStyle','Docked');
axes1 = axes('Parent',f1,'XTick',0:boundaries(1),'YTick',0:boundaries(2),'ZTick',0:boundaries(3),...
    'PlotBoxAspectRatio',[1 1 1],...
    'DataAspectRatio',[1 1 1.25],...
    'CameraViewAngle',10.339584907202,...
    'FontSize',18);
xlim(axes1,[1 boundaries(1)]-1);
ylim(axes1,[1 boundaries(2)]-1);
zlim(axes1,[1 boundaries(3)]-1);
xlabel('$v_2$','FontSize',22,'Interpreter','Latex');
ylabel('$v_3$','FontSize',22,'Interpreter','Latex');
zlabel('$v_4$','FontSize',22,'Interpreter','Latex');
view(axes1,[322.5 30]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');

switch L
    case 6
        dx = [2.4,-1,0];
        dy = [-1,3.4,0];
    case 4
        dx = [1.4,-.5,0];
        dy = [-.5,1.4,0];
    case 3
        dx = [.5,-.2,0];
        dy = [-.2,1.4,0];
end
axes1.XLabel.Position = dx;
axes1.YLabel.Position = dy;


end