function [tout] = mngdof_monge(d)
% Compute Monge coordinates at (x0,0) of surface defined by h0 below
%
% Input
%   d: degree of non-spherical terms in h0
%
% Output
%   tind: leading order terms of cubic Monge coefficients

% Symbolic variables
syms x y z x0 y0 z0 xu real
syms a [1 5] real

% Function h0: surface is level set of h0(x0,y0,z0)
h0 = z0 - (1 - y0^2 - x0^2)^(1/2) - a1*x0^d - a3*x0^(d-2)*y0^2 - a5*x0^(d-4)*y0^4

% Derivatives of h0 with respect to x0 and y0
hx0 = diff(h0,x0);

hy0 = diff(h0,y0);

hx00 = subs(hx0,[x0,y0],[xu,0]);

hy00 = subs(hy0,[x0,y0],[xu,0]);

%Normal vector

dh0 = [hx00;hy00;1];

nml = dh0/((1+hx00^2)^(1/2));
    
% Orthonormal basis P containing normal

p2 = [0;1;0];
p1 = cross(p2,nml);
P = [p1,p2,nml];

% Height of srface at (xu,0)
zu = subs(h0,[x0,y0,z0],[xu,0,0]);

% Translate coordinates to move origin to [xu,0,zu]
w = [xu,0,-zu];

% Rotate coordinates so that normal is in vertical direction
ws = [x,y,z]*P.' + w;

% Express h0 in rotated coordinates
hs = subs(h0,[x0,y0,z0],ws);

% Taylor polynomial of degree 3 for hs
ht = taylor(hs,[x,y,z],'Order',4);

% Coefficients of Taylor polynomial
[ch,th] = coeffs(ht,[x,y,z]);

%[th.',ch.']
% Extract coefficients for coordinate change to Monge form
hxxx = expand(ch(1));
hxxz = expand(ch(2));
hxx  = expand(ch(3));
hxyy = expand(ch(4));
hxzz = expand(ch(5));
hxz  = expand(ch(6));
hx   = expand(ch(7));
hyyz = expand(ch(8));
hyy  = expand(ch(9));
hzzz = expand(ch(10));
hzz  = expand(ch(11));
hz   = expand(ch(12));

% Coefficients of Monge form from h_remove_z.m
cmxxx = simplify((hxxx - (hxx*hxz)/hz)/hz);
cmxx  = simplify(hxx/hz);
cmxyy = simplify((hxyy - (hxz*hyy)/hz)/hz);
cmyyy = 0;
cmyy  = simplify(hyy/hz);
cmz   = 1;
cmxxy = 0;

% Berry-Hannay formula for index: >0 index 1/2; <0 index -1/2
% mng_ind = 3*cmxxx*cmxyy - cmxyy^2 + 3*cmxxy*cmyyy - cmyyy^2;

% Find a3(xu)
a3s = solve(cmxx - cmyy == 0,a3);

% Cubic terms of Monge form
umf = cmxx*x^2 + cmyy*y^2 + cmxxx*x^3 + cmxyy*x*y^2;

% Evaluate cubic terms of Monge form with a3 = a3s
umfs = subs(umf,a3,a3s);

% Taylor series of umfst with respect to xu and order d-2
umfst = taylor(umfs,xu,'Order',d-2)

% Coefficients of umfst as polynomial of [x,y]
[cumfst,tumfst] = coeffs(umfst,[x,y])

% Factor Berry-Hannay formula as (3*cmxxx-cmxyy)cmxyy and compute first factor
cmds = subs(3*cmxxx-cmxyy,a3,a3s);

% Taylor series of cmds and its coefficients
t3d = expand(taylor(cmds,xu,'Order',d));
[ct3d,tt3d] = coeffs(expand(t3d/a1),xu);

% Evaluate second factor cmxyy at umbilic and compute its Taylor series with coefficients
cmxyys = subs(cmxyy,a3,a3s);
txyy = expand(taylor(cmxyys,xu,'Order',d-2));
[ctxyy,ttxyy] = coeffs(expand(txyy/a1),xu);

% Compute Taylor series of the product of the two factors

tind = expand(taylor(t3d*txyy,xu,'Order',2*d-1));

% tind has factor a1^2. Compute coefficients of quotient

[ctind,ttind] = coeffs(expand(tind/a1^2),xu);

% Output leading order coefficients of factors and product

tout = [ct3d(1),ctxyy(1),ctind(1),cumfst(1)/(a1*xu^(d-3))];
