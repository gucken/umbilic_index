function [wup,wun] = mng4_umbnewtony(p) 
% Umbilic points of mng3 for parameter p
%p = [-0.001,0.001,0,0,-0.001,0]
lout = mng4_umbnewton([0,0.06],p)
wup = lout(1:2)
%wup = 0.002000730780419   0.060332712369760
lout = mng4_umbnewton([0,-0.06],p)
wun = lout(1:2)
%wun = 0.001999214309117  -0.066296524252685

