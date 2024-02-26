%Symbolic variables

syms u v w x y z hx hy hz hxx hxy hyy hxz hyz hzz a0 a1 a2 a3 b0 b1 b2 b3 b4 ld lo real

%Surface normal

Nn = [-hx;-hy;1]

%Hessian of h

d2h = [-hxx -hxy 0; -hxy -hyy 0; 0 0 0]

%Tangent space variables

vv = [u;v;w]

%Hessian image of tangent vector

d2hvv = d2h*vv

%Modified Rodrigues matrix

rodm = [Nn,vv,d2hvv]

%Modified Rodrigues matrix determinant: zeros give principal vectors

rod = det(rodm)

%w coordinate for tangent vector

ws = hx*u+hy*v

%Reexpress Rodrigues determinant

rods = expand(subs(rod,w,ws))

%Polynomial of degree 4

ha =  x^2/2+y^2/2 + ld*(x^2 -y^2) + lo*2*x*y + ...
    a3*x^3 + a2*x^2*y + a1*x*y^2 + a0*y^3 + ...
    b4*x^4 + b3*x^3*y + b2*x^2*y^2 + b1*x*y^3 + b0*y^4

%Derivatives of Monge patch

hax = diff(ha,x)

hay = diff(ha,y)

haxx = simplify(diff(hax,x))

haxy = simplify(diff(hax,y))

hayy = simplify(diff(hay,y))

%Evaluate Rodrigues criterion on Monge patch

has = simplify(subs(rods,[hx,hy,hxx,hxy,hyy],[hax,hay,haxx,haxy,hayy]))

%Taylor coefficients of Rodrigues

[ca,ta] = coeffs(has,[u,v])

cauu = ca(1)

cauv = ca(2)

cavv = ca(3)

cao = (ca(1)-ca(3))/2


