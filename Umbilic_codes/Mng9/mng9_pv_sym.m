%Symbolic variables

syms u v w x y z b1 b2 hx hy hz hxx hxy hyy hxz hyz hzz real

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

%Polynomial of degree d-1

d = 9

syms a [1 d] real
ax = poly2sym(a,x)
axy = expand(subs(y^(d-1)*ax,x,x/y))

%Monge patch for sphere + polynomial degree d

hs0 = 1 - sqrt(1 - x^2 - y^2)
ha = hs0 +axy

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