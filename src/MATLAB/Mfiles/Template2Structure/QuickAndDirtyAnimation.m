function QuickAndDirtyAnimation(P0, DP, SAbase, Bbase, knots, ...
                                scale, view_, haxis, titlestr)
% animates dynamics of interconnected structure
% INPUTS:
% P0 = control points of starting configuration of conductor (3xN matrix)
% DP = change in position of control points (3*N)x(nsteps) matrix, where
%      the control point coordinates are stored in column major form in
%      each column
% SAbase = coordinates of base of surge arrested (3x1 vector)
% Bbase = coordinates of base of bushing
% knots = knot vector
% scale = scale factor for displacements; what is plotted is P0(:)+scale*DP
% view_ = 3D graph viewpoint spec view(view_)
% haxis = handle of the axis to plot on
% titlestr = string to use for title

deformedP = P0(:) + scale*DP;

% establish axis limits
maxX = max(max(deformedP(1:3:end,:)));
minX = min(min(deformedP(1:3:end,:)));
maxY = max(max(deformedP(2:3:end,:)));
minY = min(min(deformedP(2:3:end,:)));
maxZ = max(max(deformedP(3:3:end,:)));
minZ = min(min(deformedP(3:3:end,:)));
maxX = max([maxX SAbase(1) Bbase(1)]);
minX = min([minX SAbase(1) Bbase(1)]);
maxY = max([maxY SAbase(2) Bbase(2)]);
minY = min([minY SAbase(2) Bbase(2)]);
maxZ = max([maxZ SAbase(3) Bbase(1)]);
minZ = min([minZ SAbase(3) Bbase(1)]);

xlim = [minX-0.1*(maxX-minX) maxX+0.1*(maxX-minX)];
ylim = [minY-0.1*(maxY-minY) maxY+0.1*(maxY-minY)];
zlim = [minZ-0.1*(maxZ-minZ) maxZ+0.1*(maxZ-minZ)];

% animate
N = size(deformedP,1)/3; % number of control points
for k = 1:size(deformedP,2)
    deformedP_ = reshape(deformedP(:,k),3,N);
    % plot
    bb_ = spmak(knots,deformedP_);
    points = fnplt(bb_);
    plot3(haxis,points(1,:),points(2,:),points(3,:)) % conductor
    hold(haxis,'on')
    plot3(haxis,...
          [SAbase(1) deformedP_(1,1)], ...
          [SAbase(2) deformedP_(2,1)], ...
          [SAbase(3) deformedP_(3,1)]); % deflected SA
    plot3(haxis,...
          [Bbase(1) deformedP_(1,N)], ...
          [Bbase(2) deformedP_(2,N)], ...
          [Bbase(3) deformedP_(3,N)]); % deflected bushing
    hold(haxis,'off')
    set(haxis,...
        'XLim',xlim,'YLim',ylim,'ZLim',zlim, ...
        'DataAspectRatio',[1 1 1]);
    view(haxis, view_);
    title(haxis, titlestr);
    drawnow
    %frame = getframe(gcf);
    %im{k} = frame2im(frame);
    pause(0.001)
    % fprintf('k = %d\n', k);
end

% filename = 'testAnimated.gif'; % Specify the output file name
% for idx = 1:size(deformedP,2)
%     [A,map] = rgb2ind(im{idx},256);
%     if idx == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
%     end
% end
