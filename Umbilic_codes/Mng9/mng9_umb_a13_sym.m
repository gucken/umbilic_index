
% Determine equation a3(x) for umbilics emerging from origin along x axis in family of surfaces
% h_a(x,y) = 1 - sqrt(1 - x^2 - y^2) +a1*x^8 +a2*x^7*y + a3*x^6*y^2
%
% Symbolic formula for Rodrigues criterion: ca are coefficients of R as quadratic function of (u,v)
mng9_pv_sym

% Restrict (u,v) coefficients to x axis
cas = subs(ca,y,0)

% Factor coefficients
cas1f = factor(cas(1))
cas3f = factor(cas(3))

% Coefficient of uv evaluated at [a1,a2] = [-1/28,0]
cas2 = expand(subs(cas(2),[a1,a2],[-1/28,0]))

% Factor cas2 -- x^6 is factor
cas2f = factor(cas2)

% Extract factor not divisible by x.
cas2f8 = cas2f(8)

% Evaluate at a3 = 0
a3n = subs(cas2f8,a3,0)

% Linear coefficient of a3 in cas2f8
a3d = diff(cas2f8,a3)

% Ratio -a3n/a3d approximates umbilic curve
a3x = -a3n/a3d

% Output of script
% a3x = (63*x^2*(1 - x^2)^(1/2) - 2*x^8 + 2*x^10 - 49*(1 - x^2)^(1/2))/(4*x^14*(1 - x^2)^(1/2) - 4*x^16*(1 - x^2)^(1/2) - 28*x^8 + 28*x^10 + 49*(1 - x^2)^(1/2))


