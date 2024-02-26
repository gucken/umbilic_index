p = [0,0,0,0,0,0,0,0,0,-0.01,0]
[wup,wun] = mng4_umbnewtony(p)
%
set(0, 'DefaultAxesBox','on')
fh = figure(10);
set(gcf,'Position',[100,100,864,216]);
clf
%
subplot(1,3,1)
hold on
nsteps = 1e4;
tstep = 1e-2;
stop = 0;
xlabel('x');
ylabel('y');
dp = 0.5;
[ytraj1,ye1] = mng4_traj_yc(wun + [0,1e-10],p,nsteps,tstep,stop);
plot(ytraj1(:,1),ytraj1(:,2),'b.')
[ytraj2,ye2] = mng4_traj_yc([0,0.3],p,nsteps,tstep,stop);
plot(ytraj2(:,1),ytraj2(:,2),'g.')
[ytraj3,ye3] = mng4_traj_yc([0,0.4],p,nsteps,tstep,stop);
plot(ytraj3(:,1),ytraj3(:,2),'g.')
plot(wup(:,1),wup(:,2),'k.','MarkerSize',10)
plot(wun(:,1),wun(:,2),'k.','MarkerSize',10)
axis padded;
subplot(1,3,2)
p(10) = 0;
hold on
xlabel('x');
ylabel('y');
dp = 0.5;
[ytraj4,ye4] = mng4_traj_yc([0,0.1],p,nsteps,tstep,stop);
plot(ytraj4(:,1),ytraj4(:,2),'g.')
[ytraj5,ye5] = mng4_traj_yc([0,0.3],p,nsteps,tstep,stop);
plot(ytraj5(:,1),ytraj5(:,2),'g.')
[ytraj6,ye6] = mng4_traj_yc([0,0.4],p,nsteps,tstep,stop);
plot(ytraj6(:,1),ytraj6(:,2),'g.')
plot(0,0,'k.','MarkerSize',10)
axis padded;
subplot(1,3,3)
p(10) = 0.01;
hold on
[wup,wun] = mng4_umbnewtonx(p)
xlabel('x');
ylabel('y');
dp = 0.5;
[ytraj7,ye7] = mng4_traj_yc(wun + [1e-10,0],p,nsteps,tstep,stop);
plot(ytraj7(:,1),ytraj7(:,2),'b.')
[ytraj8,ye8] = mng4_traj_yc([0.3,0],p,nsteps,tstep,stop);
plot(ytraj8(:,1),ytraj8(:,2),'g.')
[ytraj9,ye9] = mng4_traj_yc([0.4,0],p,nsteps,tstep,stop);
plot(ytraj9(:,1),ytraj9(:,2),'g.')
plot(wup(:,1),wup(:,2),'k.','MarkerSize',10)
plot(wun(:,1),wun(:,2),'k.','MarkerSize',10)
axis padded;
print(fh,  '-depsc', './umb_fig_pp0.eps')
