% Compute z = h(x,y): surface is graph of h

function [z] = mng9_z(w,p)

x = w(:,1); y = w(:,2); 
a1 = p(1); a2 = p(2); a3 = p(3); a4 = p(4); a5 = p(5); a6 = p(6); a7 = p(7); a8 = p(8); a9 = p(9);

z = 1 - (1 - y.^2 - x.^2).^(1/2) + a1.*x.^9 + a2.*x.^8.*y + a3.*x.^6.*y.^2 + a4.*x.^5.*y.^3 + a5.*x.^4.*y.^4 + a6.*x.^3.*y.^5 + a7.*x.^2.*y.^6 + a8.*x.*y.^7 + a9.*y.^8;

    
