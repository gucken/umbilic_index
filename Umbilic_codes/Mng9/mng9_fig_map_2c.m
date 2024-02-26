set(0, 'DefaultAxesBox','on')

% Values of parameters a

p = -[1/28    0   0.987143   0   0     0    2     0     1]

% Size of circle meshes

ntheta = 1000;

% Number of points in refinements of circle meshes

n0 = 100;
nthetae = ntheta + 2*n0

% Circle radii

r1 = 0.12;
r2 = 0.08;
r3 = 0.02;

%Theta mesh

thv = (2*pi*[-ntheta/2:ntheta/2-1]/ntheta)';

% Refinement meshes

n0v = linspace(-1,1,n0);
th0 = (0.01*n0v.*(abs(n0v).^3))';
thl = th0(n0/2+1:n0)-pi;
thr = th0(1:n0/2)+pi;

% Refined theta mesh 

tha = [thv;th0;thl;thr];

% Circle meshes

xyc = [cos(tha),sin(tha)];
xy1 = r1*xyc;
xy2 = r2*xyc;
xy3 = [r3*xyc(:,1)+0.1,r3*xyc(:,2)];

% Arrays that hold principal direction data along first circle 

qout1 = zeros(nthetae,3);
cdout1 = zeros(nthetae,3);
thout1 = zeros(nthetae,1);

% Compute prinicpal direction data along first circle: non-zero umbilic points are outside

for j = 1:nthetae
[cdata1,Q1,evout1] = mng9_pv_xy(xy1(j,:),p);

% Principal vector

qout1(j,:) = Q1(:,1);

% Insure angle of principal vector is in [-\pi,pi]

if evout1(1) < 0
    evout1 = -evout1;
end

% Angle of principal vector

thout1(j) = atan2(evout1(2),evout1(1));

% Quadratic coefficients in Rodrigues formula

cdout1(j,:) = cdata1;
end

% Repeat calculations along second circle: non-zero umbilic points are outside

qout2 = zeros(nthetae,3);
cdout2 = zeros(nthetae,3);
thout2 = zeros(nthetae,1);
for j = 1:nthetae
[cdata2,Q2,evout2] = mng9_pv_xy(xy2(j,:),p);
qout2(j,:) = Q2(:,1);
if evout2(1) < 0
    evout2 = -evout2;
end
thout2(j) = atan2(evout2(2),evout2(1));
cdout2(j,:) = cdata2;
end

% Repeat calculations along third circle enclosing umbilic point at (0.1,0)

qout3 = zeros(nthetae,3);
cdout3 = zeros(nthetae,3);
thout3 = zeros(nthetae,1);
for j = 1:nthetae
[cdata3,Q3,evout3] = mng9_pv_xy(xy3(j,:),p);
qout3(j,:) = Q3(:,1);
if evout3(1) < 0
    evout3 = -evout3;
end
thout3(j) = atan2(evout3(2),evout3(1));
cdout3(j,:) = cdata3;
end

% Plot angles of first principal direction as function of \theta on all three circles
% Large circle: blue
% Smaller circle around origin: green
% Circle around umbilic at (0.1,0): magenta 

fh = figure(23);
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
e = 0.01;
axis([-pi-e,pi+e,-pi/2,pi/2])
drawnow

% Plot approximate BG curves along all circles using both first and third quadratic coefficients

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

% Compute lines of curvature

% Number of lines of curvature

nthetap = 40;

% Number of time steps 

tlength = 1000;

% Time step length

tstep = 1e-3;

% Theta mesh

thp = (2*pi*[1:nthetap]/nthetap)';

% Circle mesh in (x,y) plane

xyp = r1*[cos(thp),sin(thp)];

% z values on circle mesh

zp = mng9_z(xyp,p);

% Mesh of initial points on surface

ptsin = [xyp,zp];

% Set up plot

subplot(1,3,3)
hold on
xlabel('x');
ylabel('y');
zlabel('z');
dp = 0.2;
axis([-0.2,0.2,-0.2,0.2,-0.1,0.1]);

% Plot umilic points as black dots and three circles used to compute indices

plot(0,0,'k.','MarkerSize',16,'MarkerFaceColor','k');
plot(0.1,0,'k.','MarkerSize',16,'MarkerFaceColor','k');
plot(-0.1,0,'k.','MarkerSize',16,'MarkerFaceColor','k');
plot(xy1(:,1),xy1(:,2),'b','Marker','.','MarkerSize',5,'LineStyle','none')
plot(xy2(:,1),xy2(:,2),'g','Marker','.','MarkerSize',5,'LineStyle','none')
plot(xy3(:,1),xy3(:,2),'m','Marker','.','MarkerSize',5,'LineStyle','none')


% Compute and plot lines of curvature both forward and backwards

for j = 1:nthetap
[ytraj] = mng9_traj2(xyp(j,:),p,tlength,tstep);
plot3(ytraj(:,1),ytraj(:,2),ytraj(:,3),'k.','LineWidth',0.01);
drawnow
[ytrajb] = mng9_traj2(xyp(j,:),p,tlength,-tstep);
plot3(ytrajb(:,1),ytrajb(:,2),ytrajb(:,3),'k.','LineWidth',0.01);
drawnow
end

% Print figure
print(fh,  '-depsc', './Figures_umb/mng9_fig_map2c.eps')
