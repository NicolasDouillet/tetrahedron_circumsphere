function [I, r, rc] = tetrahedron_circumsphere(A, B, C, D, option_display)
%% tetrahedron_circumsphere : function to compute, display, and save the circumsphere
% centre and radius to a given tetrahedron.
%
% Author & support : nicolas.douillet (at) free.fr, 2017-2022.
%
%
% Syntax
%
% tetrahedron_circumsphere(A, B, C, D);
% tetrahedron_circumsphere(A, B, C, D, option_display);
% [I, r, rc] = tetrahedron_circumsphere(A, B, C, D, option_display);
%
%
% Description
%
% tetrahedron_circumsphere(A, B, C, D) computes and displays the circumsphere of ABCD tetrahedron.
% tetrahedron_circumsphere(A, B, C, D, option_display) displays ABCD tetrahedron and with its circum sphere
% when option_display is set either to logical true or real numeric 1, and doesn't when it is set to logical
% false or real numeric 0.
%
% [I, r, rc] = tetrahedron_circumsphere(A, B, C, D, option_display) stores the results in [I, r, rc] vector.
%
%
% See also SPHERE
%
%
% Input arguments
%
%       [Ax]
% - A = [Ay], real column vector double, one tetrahedron vertex XYZ coordinates. Size(A) = [3,1].
%       [Az]
%
% - B, C, D : same type and description as A, here above.
%
% - option_display : logical *true(1) / false(0), to enable/disable the display mode.
%
%
% Output arguments
%
%     [Ix]
% I = [Iy], real column vector double, the circumsphere centre XYZ coordinates.
%     [Iz]
%
% - r : real scalar double, the radius of the circumsphere.
%
% - rc : logical *true(1) / false(0). The return code. rc is true when the
%        outputs are valid and false when they are invalid (degenerated cases).
%
%
% Example #1
%
% Random tetrahedron
% V = 2*(rand(3,4)-0.5);
% tetrahedron_circumsphere(V(:,1),V(:,2),V(:,3),V(:,4));
%
%
% Example #2 
%
% Regular tetrahedron in the unit sphere
% A = [2*sqrt(2)/3 0 -1/3]';
% B = [-sqrt(2)/3 sqrt(6)/3 -1/3]';
% C = [-sqrt(2)/3 -sqrt(6)/3 -1/3]';
% D = [0 0 1]';
% [I,r] = tetrahedron_circumsphere(A,B,C,D,true) % expected : I = [0 0 0]; r = 1;
%
%
% Example #3
%
% Flat / degenerated tetrahedron
% A = [0 0 0]';
% B = [1 0 0]';
% C = [0 1 0]';
% D = [1 1 0]';
% [I,r,rc] = tetrahedron_circumsphere(A,B,C,D,true);
% rc % expected : rc = 0


%% Input parsing
assert(nargin > 3, 'Not enought input arguments. Four input points are required to define one tetrahedron.');
assert(nargin < 6, 'Error : Too many input arguments.');
assert(isequal(size(A),size(B),size(C),size(D),[3 1]),'All inputs points must have the same size.');
assert(isequal(ndims(A),ndims(B),ndims(C),ndims(D),2),'All inputs points must have the same number of dimensions (2).');
assert(isreal(A) && isreal(B) && isreal(C) && isreal(D),'All inputs points must contain real numbers only.');
assert(numel(A) == 3,'Input points must have exactly 3 elements.');

if nargin < 5    
    option_display = true;            
end

%% Body

% 0 Process ABCD degenerated cases (flat or aligned points)
rc = true;

n_ABC = cross(B-A,C-A);
n_BCD = cross(C-B,D-B);

if norm(cross(n_ABC,n_BCD)) == 0    
    warning('ABCD tetrahedron is flat or degenerated; ABCD circumsphere and its centre are irrelevant.');
    rc = false;
end

V = cat(2,A,B,C,D); % tetrahedron vertex array
u12 = V(:,2)-V(:,1);
u13 = V(:,3)-V(:,1);
n123 = cross(u12,u13);
d123 = -n123(1,1)*V(1,1)-n123(2,1)*V(2,1)-n123(3,1)*V(3,1);

% Flat tetrahedron specific case
if (n123(1,1)*V(1,4)+n123(2,1)*V(2,4)+n123(3,1)*V(3,4)+d123 == 0)
    
    I = Inf*ones(3,1);
    r = Inf;
    warning('Input tetrahedron is flat -coplanar vertices-. Circumsphere centre and radius are rejected to infinity.');
    
else

% Generic tetrahedron case
    delta = det([V(:,1)' 1;...
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
    
    I = [Dx; Dy; Dz]/2/delta;
    
    r = sqrt(Dx^2+Dy^2+Dz^2-4*delta*gamma)/2/abs(delta);
end


%% Display
if option_display
        
    V = cat(2,A,B,C,D);
    [Sx,Sy,Sz] = sphere(60);
    cmap_tetra = cat(3,zeros(size(Sx)),zeros(size(Sy)),ones(size(Sz)));
    vtx_triplets = combnk(1:4,3);
    
    figure;
    plot3(I(1,1),I(2,1),I(3,1),'ko','Linewidth',4), hold on;
    s = surf(r*Sx+I(1,1),r*Sy+I(2,1),r*Sz+I(3,1),cmap_tetra); hold on; shading interp;
    alpha(s,0.4);
    
    for k = 1:size(vtx_triplets,1)
        f = fill3(V(1,vtx_triplets(k,:)),V(2,vtx_triplets(k,:)),V(3,vtx_triplets(k,:)),'r','EdgeColor','r'); hold on;
        alpha(f,0.6);
    end
    
    axis equal, axis tight;
    ax = gca;
    ax.Clipping = 'off';
    [az,el] = view(3);
    camlight(az,el);
    
end


end % tetrahedron_circumsphere