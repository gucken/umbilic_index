function [ytraj] = mng9_traj2(w,p,nsteps,del)   
 
% Step along a line of curvature 
% w: initial point
% nsteps: number of points
% del: step length

 
% Initialize storage
 
ytraj = zeros(nsteps,3);

% Limit calculations to disk of radius (0.5)^(1/2)

if w(1)^2+w(2)^2 > 0.5
    return
end
    
% Project initial point to surface 

z = mng9_z(w,p);
w = [w(1:2),z];

% Compute principal frame 
% Set orientation vector vin along maximal principal vector in right half plane

[cdata,Q,evout] = mng9_pv_xy(w,p);
vin = evout(:,1)';
if evout(1) <= 0
    evout = -evout;
end


% Step along trajectory

for j = 1:nsteps 
    
% Store current point

    ytraj(j,:) = w; 
% One step, then project back to surface

    [wout,vout] = mng9_step(@mng9_lcurv2,w,p,del,vin);  
    z = mng9_z(wout,p);
    w = [wout(1:2),z];
    if w(1)^2+w(2)^2 > 0.5
        return
    end
   
    
% Update orientation vector vin

    vin = vout;
end
