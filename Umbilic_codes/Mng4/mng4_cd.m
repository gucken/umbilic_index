function [cd] = mng4_cd(w,p)
%   Detailed explanation goes here
[hcoeff] = mng4_hcoeffs(w,p);
h = hcoeff(1);
hx = hcoeff(2);
hy = hcoeff(3);
hz = hcoeff(4);
hxx = hcoeff(5);
hxy = hcoeff(6);
hyy = hcoeff(7);


% Gradient of h and unit normal

dh = [hx;hy;hz];
nml = dh/norm(dh);
    
d2h = [hxx,hxy,0;hxy,hyy,0;0,0,0]/2;

% Orthonormal basis P containing normal

p1 = [hz;0;-hx]/sqrt(hx^2+hz^2);
p2 = cross(nml,p1);

Q = [p1,p2,nml];

% Diagonalize projected Hessan in basis Q and get cdata

d2hq = Q'*d2h*Q;

cd = [(d2hq(1,1)-d2hq(2,2))/2,d2hq(1,2)];

end