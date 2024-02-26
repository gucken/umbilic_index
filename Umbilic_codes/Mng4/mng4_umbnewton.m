function [wout] = mng4_umbnewton(w,p)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
inc = 1e-8;
for j = 1:20
    wout = [w,p];
    cd = mng4_cd(w,p)
    if cd(1)^2 + cd(2)^2 < 1e-30
        j=j
        cd = cd
        break
    end
    % Finite difference Jacobian
    dcd = [mng4_cd(w+[inc,0],p)-cd;mng4_cd(w+[0,inc],p)-cd]/inc;
    %Update point
    w = w - cd/dcd
end
