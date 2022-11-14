%% tetrahedron_circumsphere
%
% Function to compute the the circumsphere
% centre and the radius to a given tetrahedron.
%
% Author & support : nicolas.douillet (at) free.fr, 2017-2022.
%
%% Syntax
% C = tetrahedron_circumsphere(V1, V2, V3, V4);
%
% [C, r] = tetrahedron_circumsphere(V1, V2, V3, V4);
%
%% Description
% C = tetrahedron_circumsphere(V1, V2, V3, V4) computes coordinates of C, which is the
% circumsphere centre of the tetrahedron ('V1, V2, V3, V4).
%
% [C, r] = tetrahedron_circumsphere(V) also returns the circumsphere radius.
%
%% See also
%
% | <https://fr.mathworks.com/matlabcentral/fileexchange/119788-triangle-circumcircle-3d-2d?s_tid=srchtitle triangle_circumcircle> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/73490-point-to-plane-distance?s_tid=prof_contriblnk point_to_plane_distance> |
%
%% Input argument
%
%         [V1x]
% - V1 = [V1y], real column vector double, one tetrahedron vertex XYZ coordinates. Size(V1) = [3,1].
%         [V1z]
%
% - V2, V3, V4 : same type and description as V1, here above.
%
%% Output arguments
%
%      [Cx]
% C = [Cy], real column vector double, the circumsphere centre XYZ coordinates.
%      [Cz]
%
% - radius : real scalar double, the radius of the circumsphere.
%
%% Example #1
%
% Random tetrahedron
V = 2*(rand(3,4)-0.5);
[C, radius] = tetrahedron_circumsphere(V(:,1),V(:,2),V(:,3),V(:,4));
[Sx,Sy,Sz] = sphere(60);
figure;
set(gcf,'Color',[0 0 0]);
plot3(C(1,1),C(2,1),C(3,1),'ro','Linewidth',5), hold on;
couples = combnk(1:4,2);
for k = 1:size(couples,1)
    line(V(1,couples(k,:)),V(2,couples(k,:)),V(3,couples(k,:)), 'Color', [0 0 1], 'Linewidth',2), hold on;
end
surf(radius*Sx+C(1,1),radius*Sy+C(2,1),radius*Sz+C(3,1)), shading interp, 
set(gca,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1],'FontSize',16);
colormap([0 1 1]);
alpha(0.3);
camlight left;
axis equal, axis tight;

%% Example #2
%
% Flat tetrahedron
V1 = [2*sqrt(2)/3 0 -1/3]';
V2 = [-sqrt(2)/3 sqrt(6)/3 -1/3]';
V3 = [-sqrt(2)/3 -sqrt(6)/3 -1/3]';
V4 = [0 0 -1/3]';
[C, radius] = tetrahedron_circumsphere(V1,V2,V3,V4);