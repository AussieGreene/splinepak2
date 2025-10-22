%T3fcrimr Rate Version that trims the tetpart -- general PDE
% 2025: 10/9

[xp,yp,zp] = setdomain3;
ax = min(xp); bx = max(xp); ay = min(yp); by = max(yp);
az = min(zp); bz = max(zp); np = length(xp);

% choose PDE
pde = input('input pde ');
 [u,f,g,a1,a2,a3,a4,a5,a6] = setpde3(pde);

% Choose B
nb = input('input nb ');
[xb,yb,zb] = pcdown(xp,yp,zp,nb);
gb = u(xb,yb,zb);

mi = 50; %input('input m for grid of inpts ');
tic
[xr,yr,zr] = inpts3(xp,yp,zp,mi); toc;  ur = u(xr,yr,zr);


lambda = .01; lams = .01; 
d = input('input d ');  nx = input('input initial nx ');  
rmax = input('input rmax '); em = zeros(1,rmax); er = zeros(1,rmax);
h = zeros(1,rmax);


for i= 1:rmax
 ny = nx; nz = nx;  [xo,yo,zo,TETo] =  type5(nx,ny,nz,ax,bx,ay,by,az,bz);
 nto = length(TETo);
% plottet(xo,yo,zo,TETo,'b',1); nto = length(TETo);
 [x,y,z,TET,tetout] = trimtet(xo,yo,zo,TETo,xr,yr,zr); nt = length(TET);
% plottet(xo,yo,zo,TETo(tetout,:),'r',1);
 [TET,E,F,evertex,fvertex,tetvol,fstar,fb] = tetlistsb(x,y,z,TET);
% figure; plottet(x,y,z,TET,'b',1); nt = length(TET);
% tetramesh(TET,[x,y,z],'FaceAlpha',0.3,'LineWidth',2);
 h(i) = meshsize3(x,y,z,evertex);

 S = c1smooth3(d,x,y,z,TET,E,F,evertex,fvertex,fstar,fb);
 [c,G,t1,t2] = ipbc0d3(d,x,y,z,TET,E,F,evertex,fvertex,...
    a1,a2,a3,a4,a5,a6,f,g,xb,yb,zb,S,lambda,lams);
 val = valsp3(d,x,y,z,TET,E,F,evertex,fvertex,xr,yr,zr,c);
 err = u(xr,yr,zr) - val; em(i) = norm(err,inf); er(i) = erms(err);
 fprintf(' %d & %.2f & %.2f &  %d & %d & %d &  %.2e & %.2e & %.2e \\cr \n',...
     nx,t1,t2,nto,nt,length(c),condest(G), norm(err,inf),erms(err));
  nx = nx + 1;

end;

for r = 2:rmax
 fprintf('%.2f & %.2f  \\cr \n ',...
      log2(em(r-1)/em(r))/log2(h(r-1)/h(r)),...
      log2(er(r-1)/er(r))/log2(h(r-1)/h(r)));
end;

