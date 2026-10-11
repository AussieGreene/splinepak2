%T3exact % 2026: 10/7  Laplace, exact spline solution

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


d = input('input d ');
mx = input('input mx ');   my = input('input my ');   mz = input('input mz ');  
lambda = .01;   lams = .01; mc = 3; 
x = linspace(ax,bx,mx); y = linspace(ay,by,my); z = linspace(az,bz,mz);
[x,y,z,TET] =  type5(mx,my,mz,ax,bx,ay,by,az,bz);
  [TET,E,F,evertex,fvertex,tetvol,fstar,fb] = tetlistsb(x,y,z,TET);
S = c1smooth3(d,x,y,z,TET,E,F,evertex,fvertex,fstar,fb);
[c,H,t1,t2] = ipbf0d3L(d,x,y,z,TET,E,F,evertex,fvertex,f,g,xb,yb,zb,S,lambda,lams);
%[c,H,t1,t2] = ipbc0d3L(d,x,y,z,TET,E,F,evertex,fvertex,f,g,xb,yb,zb,S,lambda,lams,mc);

 val = valsp3(d,x,y,z,TET,E,F,evertex,fvertex,xr,yr,zr,c);
 err = val - u(xr,yr,zr); 
fprintf('%d & %.2f & %.2f & %d & %.2e &  %.2e &  %.2e \\cr \n',...
  mx,t1,t2,length(c),sqrt(condest(H'*H)), norm(err,inf),erms(err));
