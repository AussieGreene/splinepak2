%T3fipcimc  Test mc 
% IPBC using trivariate splines % for Laplace only
% 2025: 7/30

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

mi = 50; %input('input m for grid of inpts ');
tic
[xr,yr,zr] = inpts3(xp,yp,zp,mi); toc;  ur = u(xr,yr,zr);

d = input('input d ');  nx = input('input nx ');   lambda = .01; lams = .01;  

for mc = 2:5
 ny = nx; nz = nx;  [x,y,z,TET] =  type5(nx,ny,nz,ax,bx,ay,by,az,bz);
  [TET,E,F,evertex,fvertex,tetvol,fstar,fb] = tetlistsb(x,y,z,TET);
 S = c1smooth3(d,x,y,z,TET,E,F,evertex,fvertex,fstar,fb);
[c,G,t1,t2] = ipbc0d3L(d,x,y,z,TET,E,F,evertex,fvertex,...
    f,g,xb,yb,zb,S,lambda,lams,mc);
 val = valsp3(d,x,y,z,TET,E,F,evertex,fvertex,xr,yr,zr,c);
 err = ur - val;
io = find(isnan(err) == 1); err(io) = [];
 fprintf(' %d & %.2f & %.2f &  %d & %.2e & %.2e & %.2e \\cr \n',...
     mc,t1,t2,length(c),condest(G), norm(err,inf),erms(err));
end;

return
figure; plottet(x,y,z,TET,'b',1.4); axis equal; hold;
plot3(xp,yp,zp,'.r');
