% T3flams -- test lambda_s for IPBF with splines on tetrahedra partitions
% 2025: 10/4 
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


d = input('input d ');   nx = input('input nx ');    lambda = 1;
lams= input('input initial lams ');  
ny = nx; nz = nx;  [x,y,z,TET] =  type5(nx,ny,nz,ax,bx,ay,by,az,bz);
[TET,E,F,evertex,fvertex,tetvol,fstar,fb] = tetlistsb(x,y,z,TET);
 S = c1smooth3(d,x,y,z,TET,E,F,evertex,fvertex,fstar,fb);

for i = 1:4
 [c,G,t1,t2] = ipbf0d3L(d,x,y,z,TET,E,F,evertex,fvertex,...
    f,g,xb,yb,zb,S,lambda,lams);
 val = valsp3(d,x,y,z,TET,E,F,evertex,fvertex,xr,yr,zr,c);
 %a = find(isnan(val)); val(a) = []; xr(a) = []; yr(a) = []; zr(a) = [];
 err = u(xr,yr,zr) - val;
  fprintf(' %.3f & %.2f & %d & %.2e & %.2e & %.2e & %.2e \\cr \n',...
      lams,t1+t2,length(c),condest(G), norm(S*c),norm(err,inf),erms(err));
 lams = 10*lams; 
end;
