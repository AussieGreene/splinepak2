%TP3flam -- explore choice of lambda with IPBF method with TP splines
% 2025: 7/20

% set up domain
[xp,yp,zp] = setdomain3;
ax = min(xp); bx = max(xp); ay = min(yp); by = max(yp);
az = min(zp); bz = max(zp); np = length(xp);

% choose PDE
pde = input('input pde ');
[u,f,g] = setpde3L(pde); 
%% [u,f,g,a1,a2,a3,a4,a5,a6] = setpde3(pde);

% Choose B
nb = input('input nb ');
[xb,yb,zb] = pcdown(xp,yp,zp,nb);
gb = u(xb,yb,zb);  mc = input('input mc ');

mi = 50; %input('input mi for grid of inpts ');
tic
[xr,yr,zr] = inpts3(xp,yp,zp,mi); toc;

dx = input('input dx '); dy = dx; dz = dx;
mx = input('input mx ');  my = mx; mz= mx;

 x = linspace(ax,bx,mx); y = linspace(ay,by,my); z = linspace(az,bz,mz);
 xe = [ax*ones(1,dx),x,bx*ones(1,dx)];
 ye = [ay*ones(1,dy),y,by*ones(1,dy)];
 ze = [az*ones(1,dz),z,bz*ones(1,dz)];


lambda = input('input lambda ');
for i = 1:4
[C,G,t1,t2] = ipbctp3L(dx,dy,dz,xe,ye,ze,f,xb,yb,zb,gb,lambda,mc);
 val = valtp3v(dx,dy,dz,x,y,z,C,xr,yr,zr); nc = length(C(:));  err = val - u(xr,yr,zr); 
 fprintf('%.2f & %.2f & %d & %.2e &  %.2e &  %.2e \\cr \n',...
        lambda,t1+t2,nc,condest(G),norm(err,inf), erms(err));
lambda = 10*lambda;
end;
