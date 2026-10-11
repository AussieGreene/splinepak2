%TP3exactf % 2026: 10/7  Laplace, exact spline solution
% comment out to get ipbftp3L or ipbctp3L

% set up domain
[xp,yp,zp] = setdomain3;  % CALL WITH sphere50.dat
ax = 0; bx = 1; ay = 0; by = 1; ay = 0; by = 1; az = 0; bz = 1;

u = @(x,y,z) abs((x-.5).^3);
f = @(x,y,z) 6*abs(x-.5);
g = @(x,y,z) u(x,y,z);

% Choose B
nb = input('input nb ');
[xb,yb,zb] = pcdown(xp,yp,zp,nb);
gb = u(xb,yb,zb);

mi = 50; %input('input m for grid of inpts ');
tic
[xr,yr,zr] = inpts3(xp,yp,zp,mi); toc;


dx = input('input dx ');  dy = input('input dy ');  dz = input('input dz '); 
mx = input('input mx ');   my = input('input my ');   mz = input('input mz ');  
lambda = .01;  mc = 3; 

x = linspace(ax,bx,mx); y = linspace(ay,by,my); z = linspace(az,bz,mz);
xe = [ax*ones(1,dx),x,bx*ones(1,dx)];
ye = [ay*ones(1,dy),y,by*ones(1,dy)];
ze = [az*ones(1,dz),z,bz*ones(1,dz)];

%[C,H,t1,t2] = ipbftp3L(dx,dy,dz,xe,ye,ze,f,xb,yb,zb,gb,lambda);
[C,H,t1,t2] = ipbctp3L(dx,dy,dz,xe,ye,ze,f,xb,yb,zb,gb,lambda,mc);
 val = valtp3v(dx,dy,dz,x,y,z,C,xr,yr,zr); nc = length(C(:));
 err = val - u(xr,yr,zr); 
fprintf('%d & %.2f & %.2f & %d & %.2e &  %.2e &  %.2e \\cr \n',...
  mx,t1,t2,nc,sqrt(condest(H'*H)),norm(err,inf),erms(err));
