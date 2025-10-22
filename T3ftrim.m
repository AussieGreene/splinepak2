%T3ftrim Rate Version that trims the tetpart
% 2025: 10/6

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


lambda = .01; lams = .01; 
d = input('input d ');  nx = input('input nx ');  

 ny = nx; nz = nx;  [xo,yo,zo,TETo] =  type5(nx,ny,nz,ax,bx,ay,by,az,bz);
% plottet(xo,yo,zo,TETo,'b',1); nto = length(TETo);
 [x,y,z,TET,tetout] = trimtet(xo,yo,zo,TETo,xr,yr,zr); 
% plottet(xo,yo,zo,TETo(tetout,:),'r',1);
 [TET,E,F,evertex,fvertex,tetvol,fstar,fb] = tetlistsb(x,y,z,TET);
% figure; plottet(x,y,z,TET,'b',1); nt = length(TET);
 tetramesh(TET,[x,y,z],'FaceAlpha',0.3,'LineWidth',2);

 S = c1smooth3(d,x,y,z,TET,E,F,evertex,fvertex,fstar,fb);
 [c,G,t1,t2] = ipbf0d3L(d,x,y,z,TET,E,F,evertex,fvertex,...
    f,g,xb,yb,zb,S,lambda,lams);
 val = valsp3(d,x,y,z,TET,E,F,evertex,fvertex,xr,yr,zr,c);
 err = u(xr,yr,zr) - val;
 fprintf(' %.2f & %d & %d & %d &  %.2e & %.2e & %.2e \\cr \n',...
     t1,nto,nt,length(c),condest(G), norm(err,inf),erms(err));

