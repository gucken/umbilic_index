% This script evaluates mngdof_monge to obtain the Berry-Hannay
% expression whose sign gives the index of an umbilic at (xu,0)
% for the zero level surfaces h = 1 of
% h(x,y,z) = z - (1 - y^2 - x^2)^(1/2) - a1*x^d - a3*x^(d-2)*y^2 - a5*x^(d-4)*y^4
% at (xu,0,zu) with even degrees d in [4,20]

% Symbolic variable d
syms d

% Allocate table of data
tind_data = zeros(8,5);

% Compute the Berry-Hannay index criteria
for k = 1:8
    [tout] = mngdof_monge(2*(k+2))
    tind_data(k,:) = [2*(k+2),tout];
end

% Column vector of degrees
dv = tind_data(:,1);

% Polynomial fits to coefficients
p1 = dv + dv.^2;
p2 = -dv + 3*dv.^2/2 - dv.^3/2;
p3 = -dv.^2 + dv.^3/2 + dv.^4 -dv.^5/2;
tind_fit = [dv,p1,p2,p3]

% Confirm that fits agree with tind_data
tind_fit-tind_data(:,1:4)
