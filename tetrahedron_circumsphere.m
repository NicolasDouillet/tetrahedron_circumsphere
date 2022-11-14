function [C, radius] = tetrahedron_circumsphere(V1, V2, V3, V4)
%% tetrahedron_circumsphere : function to compute the the circumsphere
% centre and the radius to a given tetrahedron.
%
% Author & support : nicolas.douillet (at) free.fr, 2017-2022.
%
%
% Syntax
% C = tetrahedron_circumsphere(V1, V2, V3, V4);
% [C, r] = tetrahedron_circumsphere(V1, V2, V3, V4);
%
%
% Description
% C = tetrahedron_circumsphere(V1, V2, V3, V4) computes coordinates of C, which is the
% circumsphere centre of the tetrahedron (V1, V2, V3, V4).
%
% [C, r] = tetrahedron_circumsphere(V1, V2, V3, V4) also returns the circumsphere radius.
%
%
% Input arguments
%
%        [V1x]
% - V1 = [V1y], real column vector double, one tetrahedron vertex XYZ coordinates. Size(V1) = [3,1].
%        [V1z]
%
% - V2, V3, V4 : same type and description as V1, here above.
%
%
% Output arguments
%
%     [Cx]
% C = [Cy], real column vector double, the circumsphere centre XYZ coordinates.
%     [Cz]
%
% - radius : real scalar double, the radius of the circumsphere.
%
%
% Example #1
%
% Random tetrahedron
% V = 2*(rand(3,4)-0.5);
% [C, radius] = tetrahedron_circumsphere(V(:,1),V(:,2),V(:,3),V(:,4));
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
% V1 = [2*sqrt(2)/3 0 -1/3]';
% V2 = [-sqrt(2)/3 sqrt(6)/3 -1/3]';
% V3 = [-sqrt(2)/3 -sqrt(6)/3 -1/3]';
% V4 = [0 0 -1/3]';
% [C, radius] = tetrahedron_circumsphere(V1,V2,V3,V4);


%% Input parsing
assert(nargin > 0, 'Error : not enough input argument.');
assert(nargin < 5, 'Error : Too many input arguments.');
assert(isequal(size(V1),size(V2),size(V3),size(V4),[3,1]) && isreal(V1) && isreal(V2) && isreal(V3) && isreal(V4),'Inputs must be 3 x 1 real column vectors.');

%% Body
V = cat(2,V1,V2,V3,V4); % tetrahedron vertex array
u12 = V(:,2)-V(:,1);
u13 = V(:,3)-V(:,1);
n123 = cross(u12,u13);
d123 = -n123(1,1)*V(1,1)-n123(2,1)*V(2,1)-n123(3,1)*V(3,1);

% Flat tetrahedron specific case
if (n123(1,1)*V(1,4)+n123(2,1)*V(2,4)+n123(3,1)*V(3,4)+d123 == 0)
    
    C = Inf*ones(3,1);
    radius = Inf;
    warning('Input tetrahedron is flat -coplanar vertices-. Circumsphere centre and radius are rejected to infinity.');
    
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