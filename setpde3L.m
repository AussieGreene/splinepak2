% setpde3L 2025: Laplace version -- cf setpde3  for general 2nd order PDEs
% creates u,f,g for PDE 
% 2025: 7/21 -- cf setpde3kap in Scompb

function [u,f,g] = setpde3L(m);

switch m

case 1 %linear 
u = @(x,y,z) x + 2*y + 3*z; f = @(x,y,z) 0*x.^0;

case 2 % cubic
u = @(x,y,z) x.^3 + 2*y.^3 + 3*z.^3; f = @(x,y,z) 6*(x + 2*y + 3*z);

case 3 % quintic
u = @(x,y,z) x.^5 + 2*y.^5 + 3*z.^5; f = @(x,y,z) 20*(x.^3 + 2*y.^3 + 3*z.^3);

case 4 % sin5x + sin5y + sin5
u = @(x,y,z) sin(5*x) + sin(5*y) + sin(5*z); 
f = @(x,y,z) -25*(sin(5*x) + sin(5*y) + sin(5*z));

case 5 % sin5x*sin5y*sin5z
u= @(x,y,z) sin(5*x).*sin(5*y).*sin(5*z);
f = @(x,y,z) -75*u(x,y,z);

case 6 % sin(8) + sin(12) + sin(14)
u = @(x,y,z) sin(8*x) + sin(12*y) + sin(14*z);
f = @(x,y,z) -64*sin(8*x) - 144*sin(12*y) - 196*sin(14*z);

case 7  % C2 sln
u = @(x,y,z) abs((x+y-z).^3);  f = @(x,y,z) 18*abs(x+y-z);


case 8 % homog BC
u = @(x,y,z) 40*(x.^2-x).*(y.^2-y).*(z.^2-z);
f = @(x,y,z) 40*(2*((y.^2-y).*(z.^2-z) + (x.^2-x).*(z.^2-z) + (y.^2-y).*(x.^2-x)));

case 9 % degree 6
u = @(x,y,z) x.*y.^2.*z.^3;  f = @(x,y,z) 2*x.*z.^3 + 6*x.*y.^2.*z;

case 10 % 8*12*14
u = @(x,y,z) sin(8*x).*sin(12*y).* sin(14*z);
f = @(x,y,z) -404*u(x,y,z);


end;
g = @(x,y,z) u(x,y,z) ;


