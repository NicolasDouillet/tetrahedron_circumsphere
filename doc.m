%% tetrahedron_circumsphere
%
% Function to compute the the circumsphere
% centre and the radius to a given tetrahedron.
%
% Author & support : nicolas.douillet (at) free.fr, 2017-2022.
%
%% Syntax
% I = tetrahedron_circumsphere(A, B, C, D);
%
% I = tetrahedron_circumsphere(A, B, C, D, option_display);
%
% [I, r, rc] = tetrahedron_circumsphere(A, B, C, D, option_display);
%
%% Description
% tetrahedron_circumsphere(A, B, C, D) computes and displays the circumsphere of ABCD tetrahedron.
%
% tetrahedron_circumsphere(A, B, C, D, option_display) displays ABCD tetrahedron and with its circum sphere
% when option_display is set either to logical true or real numeric 1, and doesn't when it is set to logical
% false or real numeric 0.
%
% [I, r, rc] = tetrahedron_circumsphere(A, B, C, D, option_display) stores the results in [I, r, rc] vector.
%
%% See also
%
% | <https://fr.mathworks.com/matlabcentral/fileexchange/119788-triangle-circumcircle-3d-2d?s_tid=srchtitle triangle_circumcircle> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/121133-tetrahedron-insphere?s_tid=srchtitle tetrahdron_insphere> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/73490-point-to-plane-distance?s_tid=prof_contriblnk point_to_plane_distance> |
%
%% Input argument
%
%        [Ax]
% - A = [Ay], real column vector double, one tetrahedron vertex XYZ coordinates. Size(A) = [3,1].
%        [Az]
%
% - B, C, D : same type and description as A, here above.
%
% - option_display : logical *true(1) / false(0), to enable/disable the display mode.
%
%
%% Output arguments
%
%      [Ix]
% I = [Iy], real column vector double, the circumsphere centre XYZ coordinates.
%      [Iz]
%
% - r : real scalar double, the radius of the circumsphere.
%
% - rc : logical *true(1) / false(0). The return code. rc is true when the
%        outputs are valid and false when they are invalid (degenerated cases).
%
%% Example #1
%
% Random tetrahedron
V = 2*(rand(3,4)-0.5);
tetrahedron_circumsphere(V(:,1),V(:,2),V(:,3),V(:,4));

%% Example #2
%
% Regular tetrahedron in the unit sphere
A = [2*sqrt(2)/3 0 -1/3]';
B = [-sqrt(2)/3 sqrt(6)/3 -1/3]';
C = [-sqrt(2)/3 -sqrt(6)/3 -1/3]';
D = [0 0 1]';
[I,r] = tetrahedron_circumsphere(A,B,C,D,true) % expected : I = [0 0 0]; r = 1;


%% Example #3
%
% Flat / degenerated tetrahedron
A = [0 0 0]';
B = [1 0 0]';
C = [0 1 0]';
D = [1 1 0]';
[I,r,rc] = tetrahedron_circumsphere(A,B,C,D,true);
rc % expected : rc = 0