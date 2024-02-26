% Compute z = h(x,y): surface is graph of h

function [z] = mng4_z(w,p)

x = w(:,1); y = w(:,2); 
a0 = p(1); a1 = p(2); a2 = p(3); a3 = p(4);
b0 = p(5); b1 = p(6); b2 = p(7); b3 = p(8); b4 = p(9);
ld = p(10); lo = p(11);

z = x.^2/2 + y.^2/2 + ld.*(x.^2-y.^2) + lo.*2.*x.*y + ...
a3.*x.^3 + a2.*x.^2.*y + a1.*x.*y.^2 + a0.*y.^3 + ...
b4.*x.^4 + b3.*x.^3.*y + b2.*x.^2.*y.^2 + b1.*x.*y.^3 + b0.*y.^4;

    
