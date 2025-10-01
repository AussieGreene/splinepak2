%TP3ipbfr % 2025: 7/20  Rate version
% 3/13 for general 2nd order PDE % 7/20

% set up domain
[xp,yp,zp] = setdomain3;
ax = min(xp); bx = max(xp); ay = min(yp); by = max(yp);
az = min(zp); bz = max(zp); np = length(xp);

% choose PDE
pde = input('input pde ');
%%[u,f,g] = setpde3L(pde);
 [u,f,g,a1,a2,a3,a4,a5,a6] = setpde3(pde);

% Choose B
nb = input('input nb ');
[xb,yb,zb] = pcdown(xp,yp,zp,nb);
gb = u(xb,yb,zb);

mi = 50; %input('input m for grid of inpts ');
tic
[xr,yr,zr] = inpts3(xp,yp,zp,mi); toc;


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

[C,G,t1,t2] = ipbftp3(dx,dy,dz,xe,ye,ze,f,xb,yb,zb,...
    a1,a2,a3,a4,a5,a6,gb,lambda);
 val = valtp3v(dx,dy,dz,x,y,z,C,xr,yr,zr); nc = length(C(:));
 err = val - u(xr,yr,zr); em(l) = norm(err,inf); er(l) = erms(err);
 fprintf('%d & %.2f & %.2f & %d & %.2e &  %.2e &  %.2e \\cr \n',...
  mx,t1,t2,nc,condest(G),em(l),er(l));

 mx = mx+1; my = my+1; mz = mz+1;
end;

for l = 2:rmax
 fprintf('%.2f & %.2f \\cr \n ',...
      log2(em(l-1)/em(l))/log2(hr(l-1)/hr(l)),...
      log2(er(l-1)/er(l))/log2(hr(l-1)/hr(l)));
end;

return

% Error on bounding box
ngx = 51; ngy = 51; ngz = 51;
[xg,yg,zg,val] = valtp3g(dx,dy,dz,xe,ye,ze,C,ngx,ngy,ngz);
U = gridval3(xg,yg,zg,u); err =  U(:) - val(:);
fprintf('%.2e & %.2e \\cr \n',norm(err,inf), erms(err));


% plot a slice of the spline
iz = input('input level for slize iz ');
ng = ngx*ngy;   io = (iz-1)*ng;
gval = reshape(val(io+1:io+ng),ngx,ngy);
figure; surfl(xg,yg,gval'); colormap(copper);

% plot a slice of the function at the same level
%uv= reshape(U(io+1:io+ng),ngx,ngy);
%figure; surfl(xg,yg,uv'); colormap(copper);

figure;  contour(xg,yg,gval',40); hold;
%a = find(abs(zp - zg(iz)) < .01);
%plot(xp(a),yp(a),'.r')



