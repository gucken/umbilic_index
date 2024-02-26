p = [-0.002,0.002,0.001,-0.001,0,0,0,0,0,-0.01,0]
[wup,wun] = mng4_umbnewtony(p)
%
set(0, 'DefaultAxesBox','on')
fh = figure(11);
set(gcf,'Position',[100,100,864,216]);
clf
%
subplot(1,3,1)
hold on
nsteps = 1e5;
tstep = 1e-2;
stop = 0;
xlabel('x');
ylabel('y');
dp = 0.5;
[ytrajp,yep] = mng4_traj_yc(wup - [0,1e-10],p,nsteps,tstep,stop);
plot(ytrajp(:,1),ytrajp(:,2),'b.','Linewidth',0.1)
[ytrajn,yen] = mng4_traj_yc(wun + [0,1e-10],p,nsteps,tstep,stop);
plot(ytrajn(:,1),ytrajn(:,2),'g.','Linewidth',0.1)
plot(wup(:,1),wup(:,2),'k.','MarkerSize',10)
plot(wun(:,1),wun(:,2),'k.','MarkerSize',10)
axis padded;

subplot(1,3,2)
hold on
[ytraj0,ye0] = mng4_traj_yc([0.1,0],p,1.15e5,tstep,stop);
plot(ye0(:,2),'.','MarkerSize',10)
xlabel("return number")
ylabel("y")

subplot(1,3,3)
hold on
xlabel('x');
ylabel('y');
dp = 0.5;
p(10) = 0.01
[wur,wul] = mng4_umbnewtonx(p)
[ytrajr,yer] = mng4_traj_yc(wur - [1e-10,0],p,nsteps,-tstep,stop);
plot(ytrajr(:,1),ytrajr(:,2),'b.','Linewidth',0.1)
[ytrajl,yel] = mng4_traj_yc(wul + [1e-10,0],p,nsteps,tstep,stop);
plot(ytrajl(:,1),ytrajl(:,2),'g.','Linewidth',0.1)
plot(wur(:,1),wur(:,2),'k.','MarkerSize',10)
plot(wul(:,1),wul(:,2),'k.','MarkerSize',10)
axis padded;


print(fh,  '-depsc', './umb_fig_pp1.eps')
