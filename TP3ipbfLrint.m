%TP3ipbfLrint % 10/9/2026
% Edit for ForearmLink test with different interior pts

% set up domain
fname='ForearmLink.stl';
stlData = stlread(fname);
points = stlData.Points;
xo = points(:,1); yo = points(:,2); zo = points(:,3);
[xs ys zs] = scalecloud(xo,yo,zo);
TRI = stlData.ConnectivityList; nt = length(TRI);

% create point cloud
np = 50000;
pc=mesh2pc([xs,ys,zs],TRI,NumPoints=np);
pts=pc.Location;
xp = pts(:,1); yp = pts(:,2); zp = pts(:,3);

ax = min(xp); bx = max(xp); ay = min(yp); by = max(yp);
az = min(zp); bz = max(zp); np = length(xp);

% choose PDE
pde = input('input pde ');
[u,f,g] = setpde3L(pde);

% Choose B
nb = input('input nb ');
[xb,yb,zb] = pcdown(xp,yp,zp,nb);
gb = u(xb,yb,zb);

mi = 50;
[xi,yi,zi] = inpts3(xp,yp,zp,mi);
inpts=setdiff([xi,yi,zi],[xp,yp,zp],'rows','stable');
xr1=inpts(:,1); yr1=inpts(:,2); zr1=inpts(:,3);

in=inpolyhedron(TRI,[xs,ys,zs],[xr1,yr1,zr1]);
xr2=xr1(in); yr2=yr1(in); zr2=zr1(in);

figure;scatter3(xr1(in),yr1(in),zr1(in),10,'.b');axis equal; hold on;
scatter3(xr1(~in),yr1(~in),zr1(~in),10,'.r');

xr1=[xr1;xp]; yr1=[yr1;yp]; zr1=[zr1;zp];
xr2=[xr2;xp]; yr2=[yr2;yp]; zr2=[zr2;zp];

dx = input('input dx '); dy = dx; dz = dx;
mx = input('input mx ');  my = mx; mz= mx;
rmax = input('input rmax '); em = zeros(1,rmax); er = zeros(1,rmax);
hr = zeros(1,rmax);
lambda = .01; 

for l = 1:rmax
 to = cputime;
 x = linspace(ax,bx,mx); y = linspace(ay,by,my); z = linspace(az,bz,mz);
 xe = [ax*ones(1,dx),x,bx*ones(1,dx)];
 ye = [ay*ones(1,dy),y,by*ones(1,dy)];
 ze = [az*ones(1,dz),z,bz*ones(1,dz)];
 hr(l)  = x(2) - x(1);

[C,H,t1,t2] = ipbftp3L(dx,dy,dz,xe,ye,ze,f,xb,yb,zb,gb,lambda);
 val = valtp3v(dx,dy,dz,x,y,z,C,xr1,yr1,zr1); nc = length(C(:));
 err1 = val - u(xr1,yr1,zr1); em(l) = norm(err1,inf); er(l) = erms(err1);
 fprintf('%d & %.2f & %.2f & %d & %.2e &  %.2e &  %.2e \\cr \n',...
  mx,t1,t2,nc,sqrt(condest(H'*H)),em(l),er(l));

 mx = mx+1; my = my+1; mz = mz+1;
end;

for l = 2:rmax
 fprintf('%.2f & %.2f \\cr \n ',...
      log2(em(l-1)/em(l))/log2(hr(l-1)/hr(l)),...
      log2(er(l-1)/er(l))/log2(hr(l-1)/hr(l)));
end;

mx = mx-rmax; my=mx; mz=mx;
for l = 1:rmax
    to = cputime;
    x = linspace(ax,bx,mx); y = linspace(ay,by,my); z = linspace(az,bz,mz);
    xe = [ax*ones(1,dx),x,bx*ones(1,dx)];
    ye = [ay*ones(1,dy),y,by*ones(1,dy)];
    ze = [az*ones(1,dz),z,bz*ones(1,dz)];
    hr(l)  = x(2) - x(1);

    [C,H,t1,t2] = ipbftp3L(dx,dy,dz,xe,ye,ze,f,xb,yb,zb,gb,lambda);
    val = valtp3v(dx,dy,dz,x,y,z,C,xr2,yr2,zr2); nc = length(C(:));
    err1 = val - u(xr2,yr2,zr2); em(l) = norm(err1,inf); er(l) = erms(err1);
    fprintf('%d & %.2f & %.2f & %d & %.2e &  %.2e &  %.2e \\cr \n',...
        mx,t1,t2,nc,sqrt(condest(H'*H)),em(l),er(l));

    mx = mx+1; my = my+1; mz = mz+1;
end;

for l = 2:rmax
    fprintf('%.2f & %.2f \\cr \n ',...
        log2(em(l-1)/em(l))/log2(hr(l-1)/hr(l)),...
        log2(er(l-1)/er(l))/log2(hr(l-1)/hr(l)));
end;