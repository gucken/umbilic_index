p = -[0  0.145  0  1   1     0    -2     0     0.2 ]

ntheta = 2000;
n0 = 100;
r1 = 0.1;
thv = (2*pi*[-ntheta/2:ntheta/2-1]/ntheta)';
n0v = linspace(-1,1,n0);
th0 = (0.02*n0v.*(abs(n0v).^3))';
thl = th0(n0/2+1:n0)-pi;
thr = th0(1:n0/2)+pi;
tha = [thv;th0;thl;thr];
tha3 = [thv;thl;thr;th0];
xyc = [cos(tha),sin(tha)];
xy1 = r1*xyc;
%
qout1 = zeros(ntheta,3);
cdout1 = zeros(ntheta,3);
thout1 = zeros(ntheta,1);
for j = 1:ntheta+2*n0
[cdata1,Q1,evout1] = mng9_pv_xy(xy1(j,:),p);
qout1(j,:) = Q1(:,1);
cdout1(j,:) = cdata1;
if evout1(1) < 0
    evout1 = -evout1;
end
thout1(j) = atan2(evout1(2),evout1(1));
end
%
%
fh = figure(31);
set(gcf,'Position',[100,100,864,216]);
clf
%
subplot(1,3,1)
hold on
xlabel('\theta')
ylabel('\psi')
plot(tha,thout1,'b','LineStyle','none','Marker','.','MarkerSize',5)
plot([-pi,pi],[-pi,pi],'k')
plot([-pi,0],[0,pi],'k')
plot([0,pi],[-pi,0],'k')
e = 0.1;
axis([-pi-e,pi+e,-pi/2,pi/2])
drawnow
%
subplot(1,3,2)
hold on
xlabel('c_{uv}')
ylabel('(c_{vv}-c_{uu})/2')
plot((cdout1(1:ntheta,3)-cdout1(1:ntheta,1))/2,cdout1(1:ntheta,2)/2,'b','LineStyle','none','Marker','.','MarkerSize',5 )
%plot(-cdout1(1:ntheta,1),cdout1(1:ntheta,2)/2,'b','LineStyle','none','Marker','.','MarkerSize',5 )
plot(0,0,'k.','MarkerSize',16);
%axis([-4e-6,4e-6,-4e-6,4e-6])
%axis([-1.06e-6,-0.99e-6,-1.5e-8,4e-8])
drawnow
%
nthetap = 40;
tlength = 1000;
tstep = 1e-3;
thp = (2*pi*[1:nthetap]/nthetap)';
rc = 0.2;
xyp = rc*[cos(thp),sin(thp)];
zp = mng9_z(xyp,p);
ptsin = [xyp,zp];
subplot(1,3,3)
hold on
xlabel('x');
ylabel('y');
zlabel('z');
dp = 0.5;
axis([-dp,dp,-dp,dp]);
plot(ptsin(:,1),ptsin(:,2),'b.','MarkerSize',15);
%
% Compute a line of curvature, then continue a second segment
for j = 1:nthetap
[ytraj] = mng9_traj2(xyp(j,:),p,tlength,tstep);
plot3(ytraj(:,1),ytraj(:,2),ytraj(:,3),'k.','LineWidth',0.01);
drawnow
[ytrajb] = mng9_traj2(xyp(j,:),p,tlength,-tstep);
plot3(ytrajb(:,1),ytrajb(:,2),ytrajb(:,3),'k.','LineWidth',0.01);
drawnow
end
print(fh,  '-depsc', './Figures_umb/mng9_fig_map3a.eps')
