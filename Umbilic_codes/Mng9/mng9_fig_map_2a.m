set(0, 'DefaultAxesBox','on')

p = -[1/28    0   1   0   0     0    2     0     1]

ntheta = 1000;
n0 = 100;
r1 = 0.12;
r2 = 0.08;
r3 = 0.02;
thv = (2*pi*[-ntheta/2:ntheta/2-1]/ntheta)';
n0v = linspace(-1,1,n0);
th0 = (0.01*n0v.*(abs(n0v).^3))';
thl = th0(n0/2+1:n0)-pi;
thr = th0(1:n0/2)+pi;
tha = [thv;th0;thl;thr];
xyc = [cos(tha),sin(tha)];
xy1 = r1*xyc;
xy2 = r2*xyc;
xy3 = [r3*xyc(:,1)+0.1,r3*xyc(:,2)];
%
qout1 = zeros(ntheta,3);
cdout1 = zeros(ntheta,3);
thout1 = zeros(ntheta,1);
for j = 1:ntheta+2*n0
[cdata1,Q1,evout1] = mng9_pv_xy(xy1(j,:),p);
qout1(j,:) = Q1(:,1);
if evout1(1) < 0
    evout1 = -evout1;
end
thout1(j) = atan2(evout1(2),evout1(1));
cdout1(j,:) = cdata1;
end
%
qout2 = zeros(ntheta,3);
cdout2 = zeros(ntheta,3);
thout2 = zeros(ntheta,1);
for j = 1:ntheta+2*n0
[cdata2,Q2,evout2] = mng9_pv_xy(xy2(j,:),p);
qout2(j,:) = Q2(:,1);
if evout2(1) < 0
    evout2 = -evout2;
end
thout2(j) = atan2(evout2(2),evout2(1));
cdout2(j,:) = cdata2;
end
%
qout3 = zeros(ntheta,3);
cdout3 = zeros(ntheta,3);
thout3 = zeros(ntheta,1);
for j = 1:ntheta+2*n0
[cdata3,Q3,evout3] = mng9_pv_xy(xy3(j,:),p);
qout3(j,:) = Q3(:,1);
if evout3(1) < 0
    evout3 = -evout3;
end
thout3(j) = atan2(evout3(2),evout3(1));
cdout3(j,:) = cdata3;
end
%
fh = figure(21);
set(gcf,'Position',[100,100,864,216]);
clf
%
subplot(1,3,1)
hold on
xlabel('\theta')
ylabel('\psi')
plot(tha,thout1,'b','LineStyle','none','Marker','.','MarkerSize',20)
plot(tha,thout2,'g','LineStyle','none','Marker','.','MarkerSize',10)
plot(tha,thout3,'m','LineStyle','none','Marker','.','MarkerSize',5)
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
plot((cdout2(1:ntheta,3)-cdout2(1:ntheta,1))/2,cdout2(1:ntheta,2)/2,'g','LineStyle','none','Marker','.','MarkerSize',5 )
plot((cdout3(1:ntheta,3)-cdout3(1:ntheta,1))/2,cdout3(1:ntheta,2)/2,'m','LineStyle','none','Marker','.','MarkerSize',5 )
%plot(-cdout1(1:ntheta,1),cdout1(1:ntheta,2)/2,'b','LineStyle','none','Marker','.','MarkerSize',5 )
plot(0,0,'k.','MarkerSize',16);
%axis([-2e-6,2e-6,-1e-7,5e-8])
drawnow
%
nthetap = 40;
tlength = 1000;
tstep = 1e-3;
thp = (2*pi*[1:nthetap]/nthetap)';
xyp = r1*[cos(thp),sin(thp)];
zp = mng9_z(xyp,p);
ptsin = [xyp,zp];
subplot(1,3,3)
hold on
xlabel('x');
ylabel('y');
zlabel('z');
dp = 0.2;
axis([-0.2,0.2,-0.2,0.2,-0.1,0.1]);
%plot(ptsin(:,1),ptsin(:,2),'r.','MarkerSize',15);
plot(0,0,'ko','MarkerSize',5,'MarkerFaceColor','k');
% plot(0.1,0,'ko','MarkerSize',5,'MarkerFaceColor','k');
% plot(-0.1,0,'ko','MarkerSize',5,'MarkerFaceColor','k');
plot(xy1(:,1),xy1(:,2),'b','Marker','.','MarkerSize',5,'LineStyle','none')
plot(xy2(:,1),xy2(:,2),'g','Marker','.','MarkerSize',5,'LineStyle','none')
plot(xy3(:,1),xy3(:,2),'m','Marker','.','MarkerSize',5,'LineStyle','none')

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
print(fh,  '-depsc', './Figures_umb/mng9_fig_map2a.eps')
