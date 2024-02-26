function [cdata,Q,evout] = mng9_pv_xy(w,p)

% Compute principal directions of surface z = h_a(x,y)

% Input
%   (x,y): surface point is (x,y,h_a(x,y))
%   p: the parameters a

% Output
%   cdata: coefficients of Hessian d^2h
%   Q: orthonormal basis aligned with principal directions
%   evout: first principal vector

x = w(1); y = w(2); 
a1 = p(1); a2 = p(2); a3 = p(3); a4 = p(4); a5 = p(5); a6 = p(6); a7 = p(7); a8 = p(8); a9 = p(9);

% gradient of h

hax = [1, 0, x/(- x^2 - y^2 + 1)^0.5 + 8*a1*x^7 + a8*y^7 + 7*a2*x^6*y + ...
    2*a7*x*y^6 + 6*a3*x^5*y^2 + 5*a4*x^4*y^3 + 4*a5*x^3*y^4 + 3*a6*x^2*y^5];
 
hay = [0, 1, y/(- x^2 - y^2 + 1)^0.5 + a2*x^7 + 8*a9*y^7 + 2*a3*x^6*y + ...
    7*a8*x*y^6 + 3*a4*x^5*y^2 + 4*a5*x^4*y^3 + 5*a6*x^3*y^4 + 6*a7*x^2*y^5];
 
% Coefficients of d^2h as polynomial in tangent space coordinates (u,v)

cauu = (1/(1 - y^2 - x^2)^0.5 + 56*a1*x^6 + 2*a7*y^6 + ...
    x^2/(- x^2 - y^2 + 1)^1.5 + 42*a2*x^5*y + 6*a6*x*y^5 + ...
    30*a3*x^4*y^2 + 20*a4*x^3*y^3 + 12*a5*x^2*y^4)*(x/(- x^2 - y^2 + 1)^0.5 + ...
    8*a1*x^7 + a8*y^7 + 7*a2*x^6*y + 2*a7*x*y^6 + 6*a3*x^5*y^2 + ...
    5*a4*x^4*y^3 + 4*a5*x^3*y^4 + 3*a6*x^2*y^5)*(y/(- x^2 - y^2 + 1)^0.5 + ...
    a2*x^7 + 8*a9*y^7 + 2*a3*x^6*y + 7*a8*x*y^6 + 3*a4*x^5*y^2 + ...
    4*a5*x^4*y^3 + 5*a6*x^3*y^4 + 6*a7*x^2*y^5) - 7*a2*x^6 - 7*a8*y^6 - ...
    (7*a2*x^6 + 7*a8*y^6 + 12*a3*x^5*y + 12*a7*x*y^5 + 15*a4*x^4*y^2 + ...
    16*a5*x^3*y^3 + 15*a6*x^2*y^4 + (x*y)/(- x^2 - y^2 + 1)^1.5)*(x/(- x^2 ...
    - y^2 + 1)^0.5 + 8*a1*x^7 + a8*y^7 + 7*a2*x^6*y + 2*a7*x*y^6 + ...
    6*a3*x^5*y^2 + 5*a4*x^4*y^3 + 4*a5*x^3*y^4 + 3*a6*x^2*y^5)^2 - ...
    12*a3*x^5*y - 12*a7*x*y^5 - 15*a4*x^4*y^2 - 16*a5*x^3*y^3 - ...
    15*a6*x^2*y^4 - (x*y)/(- x^2 - y^2 + 1)^1.5;
 
cauv = 56*a1*x^6 - 2*a3*x^6 + 2*a7*y^6 - 56*a9*y^6 + (1/(1 - y^2 - x^2)^0.5 + 56*a1*x^6 + 2*a7*y^6 + x^2/(- x^2 - y^2 + 1)^1.5 + 42*a2*x^5*y + 6*a6*x*y^5 + 30*a3*x^4*y^2 + 20*a4*x^3*y^3 + 12*a5*x^2*y^4)*(y/(- x^2 - y^2 + 1)^0.5 + a2*x^7 + 8*a9*y^7 + 2*a3*x^6*y + 7*a8*x*y^6 + 3*a4*x^5*y^2 + 4*a5*x^4*y^3 + 5*a6*x^3*y^4 + 6*a7*x^2*y^5)^2 - (1/(1 - y^2 - x^2)^0.5 + 2*a3*x^6 + 56*a9*y^6 + y^2/(- x^2 - y^2 + 1)^1.5 + 6*a4*x^5*y + 42*a8*x*y^5 + 12*a5*x^4*y^2 + 20*a6*x^3*y^3 + 30*a7*x^2*y^4)*(x/(- x^2 - y^2 + 1)^0.5 + 8*a1*x^7 + a8*y^7 + 7*a2*x^6*y + 2*a7*x*y^6 + 6*a3*x^5*y^2 + 5*a4*x^4*y^3 + 4*a5*x^3*y^4 + 3*a6*x^2*y^5)^2 + x^2/(- x^2 - y^2 + 1)^1.5 - y^2/(- x^2 - y^2 + 1)^1.5 + 42*a2*x^5*y - 6*a4*x^5*y + 6*a6*x*y^5 - 42*a8*x*y^5 + 30*a3*x^4*y^2 + 20*a4*x^3*y^3 + 12*a5*x^2*y^4 - 12*a5*x^4*y^2 - 20*a6*x^3*y^3 - 30*a7*x^2*y^4;

cavv = (7*a2*x^6 + 7*a8*y^6 + 12*a3*x^5*y + 12*a7*x*y^5 + 15*a4*x^4*y^2 + 16*a5*x^3*y^3 + 15*a6*x^2*y^4 + (x*y)/(- x^2 - y^2 + 1)^1.5)*(y/(- x^2 - y^2 + 1)^0.5 + a2*x^7 + 8*a9*y^7 + 2*a3*x^6*y + 7*a8*x*y^6 + 3*a4*x^5*y^2 + 4*a5*x^4*y^3 + 5*a6*x^3*y^4 + 6*a7*x^2*y^5)^2 + 7*a2*x^6 + 7*a8*y^6 - (1/(1 - y^2 - x^2)^0.5 + 2*a3*x^6 + 56*a9*y^6 + y^2/(- x^2 - y^2 + 1)^1.5 + 6*a4*x^5*y + 42*a8*x*y^5 + 12*a5*x^4*y^2 + 20*a6*x^3*y^3 + 30*a7*x^2*y^4)*(x/(- x^2 - y^2 + 1)^0.5 + 8*a1*x^7 + a8*y^7 + 7*a2*x^6*y + 2*a7*x*y^6 + 6*a3*x^5*y^2 + 5*a4*x^4*y^3 + 4*a5*x^3*y^4 + 3*a6*x^2*y^5)*(y/(- x^2 - y^2 + 1)^0.5 + a2*x^7 + 8*a9*y^7 + 2*a3*x^6*y + 7*a8*x*y^6 + 3*a4*x^5*y^2 + 4*a5*x^4*y^3 + 5*a6*x^3*y^4 + 6*a7*x^2*y^5) + 12*a3*x^5*y + 12*a7*x*y^5 + 15*a4*x^4*y^2 + 16*a5*x^3*y^3 + 15*a6*x^2*y^4 + (x*y)/(- x^2 - y^2 + 1)^1.5;
   
% Compute principal directions as eigenvectors of matrix A constructed from Hessian

A = [[-cauv/2,-cavv];[cauu,cauv/2]];
[V,D] = eig(A);

% Compute principal vectors 

pd1m = (V(1,1)*hax+V(2,1)*hay);
pd2m = (V(1,2)*hax+V(2,2)*hay);

% Normalize principal vectors

pv1m = pd1m/norm(pd1m);
pv2m = pd2m/norm(pd2m);

% Compute normal vector

N = cross(pv1m,pv2m);

% Orthonormal basis with first basis vector in maximal prinicipal direction

Q=[pv1m',pv2m',N'];
if D(2,2) > 0
    evout = pv2m;
else
    evout = pv1m;
end

% Coefficients of d^2h as polynomial of (u,v)

cdata=[cauu;cauv;cavv];

