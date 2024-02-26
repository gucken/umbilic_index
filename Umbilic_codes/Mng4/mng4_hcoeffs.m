function [hcoeff] = mng4_hcoeffs(w,p)

x = w(1);y=w(2);

a0 = p(1); a1 = p(2); a2 = p(3); a3 = p(4);
b0 = p(5); b1 = p(6); b2 = p(7); b3 = p(8); b4 = p(9);
ld = p(10); lo = p(11);

h = ld*(- y^2 + x^2) + a3*x^3 + a0*y^3 + b4*x^4 + b0*y^4 + x^2/2 + y^2/2 + ...
    a1*x*y^2 + a2*x^2*y + b1*x*y^3 + b3*x^3*y + b2*x^2*y^2 + 2*lo*x*y;
 
 
hx = 2*lo*y + a1*y^2 + b1*y^3 + x + 2*ld*x + 2*a2*x*y + 2*b2*x*y^2 + ...
    3*a3*x^2 + 3*b3*x^2*y + 4*b4*x^3;
 
 
hy = y - 2*ld*y + 3*a0*y^2 + 4*b0*y^3 + 2*lo*x + 2*a1*x*y + 3*b1*x*y^2 + ...
a2*x^2 + 2*b2*x^2*y + b3*x^3;

hz = 1; 
 
hxx = 1 + 2*ld + 2*a2*y + 2*b2*y^2 + 6*a3*x + 6*b3*x*y + 12*b4*x^2;
 
 
hxy = 2*lo + 2*a1*y + 3*b1*y^2 + 2*a2*x + 4*b2*x*y + 3*b3*x^2;
 
 
hyy = 1 - 2*ld + 6*a0*y + 12*b0*y^2 + 2*a1*x + 6*b1*x*y + 2*b2*x^2;

hcoeff = [h,hx,hy,hz,hxx,hxy,hyy];

%dp_coeff = [hld,hlo,ha0,ha1,ha2,ha3,hb0,hb1,hb2,hb3,hb4];
