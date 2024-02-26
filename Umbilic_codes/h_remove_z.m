% Symbolic computation of Monge form for surface defined as level set of
% Cubic polynomial h(x,y,z) with dh(0,0,0) = (0,0,az)
% Transform the cubic terms by multiplying with quadratic function 1 + g(x,y,z)
% with g(0,0,0) = 0
% Both h and g are assumed to be even in y

%Variables
syms x y z az axx axz ayy axyy ayz azz axxx axxxz axxz axzz ayyz ayyy azzz 
syms bx bz bxx bxz bzz byy

% Function h
h = az*z + axx*x^2 +axz*x*z + azz*z^2 + axxx*x^3 + axxz*x^2*z + axzz*x*z^2 ...
    + azzz*z^3 + axyy*x*y^2 + ayyy*y^3 + ayyz*y^2*z +ayy*y^2

% Function g 
g = bx*x + bz*z + bxx*x^2 + bxz*x*z + bzz*z^2 + byy*y^2 

% Taylor series of product gh
hb = taylor((1+g)*h,[x,y,z],'Order',4)

% Coefficients of Taylor series
[ch,th] = coeffs(hb,[x,y,z])

% Transpose array of coefficients 
[th.',ch.']

% Output
% [  x^3,                   axxx + axx*bx]
% [x^2*z, az*bxx + axz*bx + axxz + axx*bz]
% [  x^2,                             axx]
% [x*y^2,                   ayy*bx + axyy]
% [x*z^2, azz*bx + az*bxz + axzz + axz*bz]
% [  x*z,                     az*bx + axz]
% [  y^3,                            ayyy]
% [y^2*z,          az*byy + ayyz + ayy*bz]
% [  y^2,                             ayy]
% [  z^3,          azzz + azz*bz + az*bzz]
% [  z^2,                     azz + az*bz]
% [    z,                              az]

% Solve for [bx,bz,bxx,bxz,bzz,byy] to kill coefficients that do not appear
% in Monge form and then make substitutions in (1+g)h
gs = expand(subs(g,[bx,bz,bxx,bxz,bzz,byy],[-axz/az,-azz/az,- axxz/az + ...
    axz*axz/az^2 + axx*azz/az^2,- axzz/az + azz*axz/az^2 + axz*azz/az^2,...
    - azzz/az + azz*azz/az^2, -ayyz/az + ayy*azz/az^2]))

% Taylor expansion of (1+gs)h
hbs =expand((1+gs)*h);
hbst = taylor(hbs,[x,y,z],'Order',4)

% Coefficients of hbst
[chbs,thbs] = coeffs(hbst,[x,y,z])

% Transpose array of coefficients
[thbs.', chbs.']

% Output gives the coefficients of Monge form
% [  x^3, axxx - (axx*axz)/az]
% [  x^2,                 axx]
% [x*y^2, axyy - (axz*ayy)/az]
% [  y^3,                ayyy]
% [  y^2,                 ayy]
% [    z,                  az]
