% Computaioon of curve of umbilic points emerging from origin at parameters p 

p = [-1/28    0    -1   0   0     0    -2     0     -1];

% Inline function defining a3(x) -- output of mng9_a13_sym.m

a3x = @(x)((63*x^2*(1 - x^2)^(1/2) - 2*x^8 + 2*x^10 - 49*(1 - x^2)^(1/2))/(4*x^14*(1 - x^2)^(1/2) - 4*x^16*(1 - x^2)^(1/2) - 28*x^8 + 28*x^10 + 49*(1 - x^2)^(1/2)))

% vector of values for x

xv = 0.01*[-10:10];

% Vector for a3 values

a3v = zeros(1,21);

% Compute values for a3

for k = 1:21
    a3v(k) = a3x(xv(k));
end

% Plot curve a3 vs x

figure(4)
clf 
hold on
plot(xv,a3v+1,'.')

% Compute quadratic approximation to a3(v)

cf = polyfit(xv,a3v,2)
xf = polyval(cf,xv);

% Plot appoximation

plot(xv,xf+1)

% Table of values and differences

[a3v;xf;a3v-xf;-1-a3v]'
