
function [ytraj,ye] = mng4_traj_yc(w,p,nsteps,del,stop)
 
 % Step along a line of curvature 
 % w: initial point
 % p: parameter
 % nsteps: number of points
 % del: step length
 
% Initialize storage
 
ytraj = zeros(nsteps,3);
ye = [];

% Project initial point to surface h_a=0

w(3) = mng4_z(w,p);

if w(1)^2+w(2)^2 > 1
    return
end
 
% Compute principal frame
% Set orientation vector vin along maximal principal direction

vin = [1,0,0];

% Step along trajectory

for j = 1:nsteps 
    
    % Store current point

    ytraj(j,:) = w;

   % Adapt step size based on distance to umbilic points
    [cdata,Q,evout] = mng4_pv_xy(w,p);
    cud = sqrt((cdata(1)-cdata(3))^2/4+cdata(2)^2);
   
    % One step, then project back to surface

    [wout,vout] = mng4_step(@mng4_lcurv2,w,p,del*cud/(1+cud),vin);
    z = mng4_z(wout,p);
    w = [wout(1:2),z];
    if w(1)^2+w(2)^2 > 0.5
        ytraj  = ytraj(1:j,:);
        jend = j
        return
    end
   
    % Update orientation vector vin

    vin = vout;
    
    % Test crossing
    
    if ytraj(j,1) < 0 & w(1) > 0
        ycross = (w(1)*ytraj(j,:) - ytraj(j,1)*w)/(w(1)-ytraj(j,1));
        ye = [ye;ycross];
        if stop > 0
            ytraj = [ytraj(1:j,:);ye];
            j = j
        return
        end
    end
end
