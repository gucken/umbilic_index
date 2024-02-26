function [cdata,Q,evout] = mng4_pv_xy(w,p)

% Compute principal directions of surface z = h(x,y)

% Input
%   (x,y): surface point is (x,y,h_a(x,y))
%   p: the parameters a

% Output
%   cdata: coefficients of Hessian d^2h
%   Q: orthonormal basis aligned with principal directions
%   evout: first principal vector

hcoeff = mng4_hcoeffs(w,p);
h = hcoeff(1);
hx = hcoeff(2);
hy = hcoeff(3);
hz = hcoeff(4);
hxx = hcoeff(5);
hxy = hcoeff(6);
hyy = hcoeff(7);

% Gradient of h and unit normal
dh = [hx,hy,hz];
nml = dh/norm(dh);

%Hessian
d2h = [hxx,hxy,0;hxy,hyy,0;0,0,0]/2;

% Coefficients of d^2h as polynomial in tangent space coordinates (u,v)
x = w(1); y = w(2); 
a0 = p(1); a1 = p(2); a2 = p(3); a3 = p(4); 
b0 = p(5); b1 = p(6); b2 = p(7); b3 = p(8); b4 = p(9);
ld = p(10); lo = p(11);

cauu = (1 + 2*ld + 2*a2*y + 2*b2*y^2 + 6*a3*x + 6*b3*x*y + 12*b4*x^2)*(2*lo*y + ...
    a1*y^2 + b1*y^3 + x + 2*ld*x + 2*a2*x*y + 2*b2*x*y^2 + 3*a3*x^2 + 3*b3*x^2*y + ...
    4*b4*x^3)*(y - 2*ld*y + 3*a0*y^2 + 4*b0*y^3 + 2*lo*x + 2*a1*x*y + 3*b1*x*y^2 + ...
    a2*x^2 + 2*b2*x^2*y + b3*x^3) - 2*a2*x - 2*a1*y - 3*b3*x^2 - 3*b1*y^2 - ...
    (2*lo + 2*a1*y + 3*b1*y^2 + 2*a2*x + 4*b2*x*y + 3*b3*x^2)*(2*lo*y + a1*y^2 + ...
    b1*y^3 + x + 2*ld*x + 2*a2*x*y + 2*b2*x*y^2 + 3*a3*x^2 + 3*b3*x^2*y + 4*b4*x^3)^2 - 2*lo - 4*b2*x*y;
 
 
cauv = 4*ld - 2*a1*x + 6*a3*x - 6*a0*y + 2*a2*y - 2*b2*x^2 + 12*b4*x^2 - 12*b0*y^2 + 2*b2*y^2 - ...
    (1 - 2*ld + 6*a0*y + 12*b0*y^2 + 2*a1*x + 6*b1*x*y + 2*b2*x^2)*(2*lo*y + a1*y^2 + b1*y^3 + ...
    x + 2*ld*x + 2*a2*x*y + 2*b2*x*y^2 + 3*a3*x^2 + 3*b3*x^2*y + 4*b4*x^3)^2 + (1 + 2*ld + ...
    2*a2*y + 2*b2*y^2 + 6*a3*x + 6*b3*x*y + 12*b4*x^2)*(y - 2*ld*y + 3*a0*y^2 + 4*b0*y^3 + ...
    2*lo*x + 2*a1*x*y + 3*b1*x*y^2 + a2*x^2 + 2*b2*x^2*y + b3*x^3)^2 - 6*b1*x*y + 6*b3*x*y;
 
 
cavv = 2*lo + 2*a2*x + 2*a1*y + 3*b3*x^2 + 3*b1*y^2 + (2*lo + 2*a1*y + 3*b1*y^2 + 2*a2*x + ...
    4*b2*x*y + 3*b3*x^2)*(y - 2*ld*y + 3*a0*y^2 + 4*b0*y^3 + 2*lo*x + 2*a1*x*y + 3*b1*x*y^2 + ...
    a2*x^2 + 2*b2*x^2*y + b3*x^3)^2 - (1 - 2*ld + 6*a0*y + 12*b0*y^2 + 2*a1*x + 6*b1*x*y + ...
    2*b2*x^2)*(2*lo*y + a1*y^2 + b1*y^3 + x + 2*ld*x + 2*a2*x*y + 2*b2*x*y^2 + 3*a3*x^2 + ...
    3*b3*x^2*y + 4*b4*x^3)*(y - 2*ld*y + 3*a0*y^2 + 4*b0*y^3 + 2*lo*x + 2*a1*x*y + 3*b1*x*y^2 + ...
    a2*x^2 + 2*b2*x^2*y + b3*x^3) + 4*b2*x*y;
 
% Compute principal directions as eigenvectors of matrix A constructed from Hessian

A = [[-cauv/2,-cavv];[cauu,cauv/2]];
[V,D] = eig(A);

% Compute principal vectors 

pd1m = (V(1,1)*[1;0;-hx]+V(2,1)*[0;1;-hy]);
pd2m = (V(1,2)*[1;0;-hx]+V(2,2)*[0;1;-hy]);

% Normalize principal vectors

pv1m = pd1m/norm(pd1m);
pv2m = pd2m/norm(pd2m);

% Orthonormal basis with first basis vector in maximal prinicipal direction

Q=[pv1m,pv2m,nml'];

P = Q'*d2h*Q;

if P(1,1) > P(2,2)
    evout = pv1m';
else
    evout = pv2m';
end

% Coefficients of d^2h as polynomial of (u,v)

cdata=[cauu;cauv;cavv];


