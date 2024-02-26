function [yout,vout] = mng9_step(f,y,p,del,vin)

% Runge-Kutta step of length h
% f evaluates principal directions, using vin to set orientation
% y is input point
% del is step length
% yout is output point
% vout can be used to update vin in computing lines of curvature

    s1 = f(y,p,vin);
    s2 = f(y + del * s1/2,p,vin);
    s3 = f(y + del * s2/2,p,vin);
    s4 = f(y+ del * s3,p,vin);
    yout = y + (s1 + s2+s2 + s3+s3 + s4) * del/6;
    vout = s4;
end
