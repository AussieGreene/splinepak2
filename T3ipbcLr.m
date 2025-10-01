% T3ipbcLr -- Rate version
%2025: 7/30
% set up domain
[xp,yp,zp] = setdomain3;
ax = min(xp); bx = max(xp); ay = min(yp); by = max(yp);
az = min(zp); bz = max(zp); np = length(xp);

% choose PDE
pde = input('input pde ');
[u,f,g] = setpde3L(pde);

% Choose B
nb = input('input nb ');
[xb,yb,zb] = pcdown(xp,yp,zp,nb);
gb = u(xb,yb,zb); 
mc = input('input mc ');
mi = 50; %input('input m for grid of inpts ');
tic
[xr,yr,zr] = inpts3(xp,yp,zp,mi); toc;  ur = u(xr,yr,zr);


d = input('input d ');   nx = input('input nx ');    lambda = .01; lams = .01;  
rmax = input('input rmax '); em = zeros(1,rmax); er = zeros(1,rmax);

for l = 1:rmax
 ny = nx; nz = nx;  [x,y,z,TET] =  type5(nx,ny,nz,ax,bx,ay,by,az,bz);
  [TET,E,F,evertex,fvertex,tetvol,fstar,fb] = tetlistsb(x,y,z,TET);
 [h(l),el] = meshsize3(x,y,z,evertex);
 S = c1smooth3(d,x,y,z,TET,E,F,evertex,fvertex,fstar,fb);
[c,G,t1,t2] = ipbc0d3L(d,x,y,z,TET,E,F,evertex,fvertex,...
      f,g,xb,yb,zb,S,lambda,lams,mc);
 val = valsp3(d,x,y,z,TET,E,F,evertex,fvertex,xr,yr,zr,c);
 a = find(isnan(val)); val(a) = []; xr(a) = []; yr(a) = []; zr(a) = [];
 err = u(xr,yr,zr) - val;
 ur = u(xr,yr,zr); io = find(isnan(val) == 1);
ur(io) = []; val(io) = [];  err = ur - val;
 em(l) = norm(err,inf); er(l) = erms(err);
 fprintf(' %d & %.2f & %.2f &  %d & %.2e & %.2e & %.2e \\cr \n',...
     nx,t1,t2,length(c),condest(G), norm(err,inf),erms(err));
nx = nx + 1;
end;

for r = 2:rmax
 fprintf('%.2f & %.2f  \\cr \n ',...
      log2(em(r-1)/em(r))/log2(h(r-1)/h(r)),...
      log2(er(r-1)/er(r))/log2(h(r-1)/h(r)));
end;


return
figure; plottet(x,y,z,TET,'b',1.4); axis equal; hold;
plot3(xp,yp,zp,'.r');
