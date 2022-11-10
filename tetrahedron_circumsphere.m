function [C, radius] = tetrahedron_circumsphere(V)
%% tetrahedron_circumsphere : function to compute
% the centre and the radius of the sphere circumscribed to
% a given tetrahedron.
%
% Author & support : nicolas.douillet (at) free.fr, 2017-2022.
%
%
% Syntax
% C = tetrahedron_circumsphere(V);
% [C, r] = tetrahedron_circumsphere(V);
%
%
% Description
% C = tetrahedron_circumsphere(V) computes coordinates of C, which is the
% centre of the circumscribed sphere to the tetrahedron V.
%
% [C, r] = tetrahedron_circumsphere(V) also returns the radius of the
% circumscribed sphere.
%
%
% Input argument
%
%       [V1x V2x V3x V4x]
% - V = [V1y V2y V3y V4y], real matrix double, the ctetrahedron vertex XYZ coordinates. Size(V) = [3,4].
%       [V1z V2z V3z V4z]       
%
%
% Output arguments
%
%     [Cx]
% C = [Cy], real column vector double, the circumscribed centre XYZ coordinates.
%     [Cz]
%
% - radius : real scalar double, the radius of the circumscribed sphere.
%
%
% Example #1
%
% Random tetrahedron
% V = 2*(rand(3,4)-0.5);
% [C, radius] = tetrahedron_circumsphere(V);
% [Sx,Sy,Sz] = sphere(60);
% figure;
% set(gcf,'Color',[0 0 0]);
% plot3(C(1,1),C(2,1),C(3,1),'ro','Linewidth',5), hold on;
% couples = combnk(1:4,2);
% for k = 1:size(couples,1)
%     line(V(1,couples(k,:)),V(2,couples(k,:)),V(3,couples(k,:)), 'Color', [0 1 0], 'Linewidth',2), hold on;
% end
% surf(radius*Sx+C(1,1),radius*Sy+C(2,1),radius*Sz+C(3,1)), shading interp, 
% set(gca,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1],'FontSize',16);
% colormap([0 1 1]);
% alpha(0.3);
% camlight left;
% axis equal, axis tight;
%
%
% Example #2
%
% Flat tetrahedron
% V1 = [2*sqrt(2)/3 0 -1/3];
% V2 = [-sqrt(2)/3 sqrt(6)/3 -1/3];
% V3 = [-sqrt(2)/3 -sqrt(6)/3 -1/3];
% V4 = [0 0 -1/3];
% V = [V1' V2' V3' V4'];
% [C, radius] = tetrahedron_circumsphere(V);


%% Input parsing
assert(nargin > 0, 'Error : not enough input argument.');
assert(nargin < 2, 'Error : Too many input arguments.');
assert(nargout < 3, 'Error : too many output arguments.');
assert(size(V,1) == 3 && size(V,2) == 4, 'Error : input V must be a 3 x 4 matrix of finite coordinates.');


%% Body
u12 = V(:,2)-V(:,1);
u13 = V(:,3)-V(:,1);
n123 = cross(u12,u13);
d123 = -n123(1,1)*V(1,1)-n123(2,1)*V(2,1)-n123(3,1)*V(3,1);

% Flat tetrahedron specific case
if (n123(1,1)*V(1,4)+n123(2,1)*V(2,4)+n123(3,1)*V(3,4)+d123 == 0)
    C = Inf*ones(3,1);
    radius = Inf;
    warning('Input tetrahedron is flat -coplanar vertices-. Circumscribed centre and radius are rejected to infinity.');
else

% Generic tetrahedron case
    alpha = det([V(:,1)' 1;...
                 V(:,2)' 1;...
                 V(:,3)' 1;...
                 V(:,4)' 1]);
    
    gamma = det([sum((V(:,1)').^2,2) V(1,1) V(2,1) V(3,1);...
                 sum((V(:,2)').^2,2) V(1,2) V(2,2) V(3,2);...
                 sum((V(:,3)').^2,2) V(1,3) V(2,3) V(3,3);...
                 sum((V(:,4)').^2,2) V(1,4) V(2,4) V(3,4)]);
    
    Dx = det([sum((V(:,1)').^2,2) V(2,1) V(3,1) 1;...
              sum((V(:,2)').^2,2) V(2,2) V(3,2) 1;...
              sum((V(:,3)').^2,2) V(2,3) V(3,3) 1;...
              sum((V(:,4)').^2,2) V(2,4) V(3,4) 1]);
    
    Dy = -det([sum((V(:,1)').^2,2) V(1,1) V(3,1) 1;...
               sum((V(:,2)').^2,2) V(1,2) V(3,2) 1;...
               sum((V(:,3)').^2,2) V(1,3) V(3,3) 1;...
               sum((V(:,4)').^2,2) V(1,4) V(3,4) 1]);
    
    Dz = det([sum((V(:,1)').^2,2) V(1,1) V(2,1) 1;...
              sum((V(:,2)').^2,2) V(1,2) V(2,2) 1;...
              sum((V(:,3)').^2,2) V(1,3) V(2,3) 1;...
              sum((V(:,4)').^2,2) V(1,4) V(2,4) 1]);
    
    C = [Dx; Dy; Dz]/2/alpha;
    
    radius = sqrt(Dx^2+Dy^2+Dz^2-4*alpha*gamma)/2/abs(alpha);
end


end % tetrahedron_circumsphere