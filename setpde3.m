
% setpde3 2025: 4/27
% define coefs for a PDE and sets up a BVP Lu= f
% where Lu = a1*uxx + a2*uyy + a3*uzz + 2*a4*uxy + 2*a5*uxz + 2*a6*uyz 

function [u,f,g,a1,a2,a3,a4,a5,a6] = setpde3(pde);



a1 = @(x,y,z) x.^0;  a2 = @(x,y,z) y.^0;  a3 = @(x,y,z) z.^0;  
a4 = @(x,y,z) 0*x.^0;  a5 = @(x,y,z) 0*y.^0;  a6 = @(x,y,z) 0*z.^0; 

switch pde

case 1 %linear 
u = @(x,y,z) x + 2*y + 3*z; 
ux = @(x,y,z) x.^0; uy = @(x,y,z) 2*y.^0; uz = @(x,y,z) 3*z^0;
uxx = @(x,y,z) 0*x.^3; uyy = @(x,y,z) 0*y.^3; uzz = @(x,y,z) 0*z.^3;
uxy = @(x,y,z) 0*x.^0; uxz = @(x,y,z) 0*x.^0; uyz = @(x,y,z) 0*x.^0;

case 2 % cubic
u = @(x,y,z) x.^3 + 2*y.^3 + 3*z.^3; 
ux = @(x,y,z) 3*x.^2; uy = @(x,y,z) 6*y.^2; uz = @(x,y,z) 9*z.^2;
uxx = @(x,y,z) 6*x.^0; uyy = @(x,y,z) 12*y.^0; uzz = @(x,y,z) 18*z.^0; 
uxy = @(x,y,z) 0*x.^0; uxz = @(x,y,z) 0*x.^0; uyz = @(x,y,z) 0*x.^0;


case 3 % quintic
u = @(x,y,z) x.^5 + 2*y.^5 + 3*z.^5; 
ux = @(x,y,z) 5*x.^4; uy = @(x,y,z) 10*y.^4; uz = @(x,y,z) 15*z.^4;
uxx = @(x,y,z) 20*x.^3; uyy = @(x,y,z) 40*y.^3; uzz = @(x,y,z) 60*z.^3; 
uxy = @(x,y,z) 0*x.^0; uxz = @(x,y,z) 0*x.^0; uyz = @(x,y,z) 0*x.^0;

case 4 % sin5x + sin5y + sin5
u = @(x,y,z) sin(5*x) + sin(5*y) + sin(5*z); 
ux = @(x,y,z) 5*cos(5*x); uy = @(x,y,z) 5*cos(5*y); uz = @(x,y,z) 5*cos(5*z);
uxx = @(x,y,z) -25*sin(5*x); uyy = @(x,y,z) -25*sin(5*y); uzz = @(x,y,z) -25*sin(5*z);
uxy = @(x,y,z) 0*x.^0; uyx = @(x,y,z) 0*x.^0;
uxz = @(x,y,z) 0*x.^0; uzx = @(x,y,z) 0*x.^0;
uyz = @(x,y,z) 0*x.^0; uzy = @(x,y,z) 0*x.^0;

case 5 % sin5x*sin5y*sin5z
u= @(x,y,z) sin(5*x).*sin(5*y).*sin(5*z);
ux = @(x,y,z) 5*cos(5*x).*sin(5*y).*sin(5*z);
uy = @(x,y,z) 5*sin(5*x).*cos(5*y).*sin(5*z);
uz = @(x,y,z) 5*sin(5*x).*sin(5*y).*cos(5*z);
uxx = @(x,y,z) -25*u(x,y,z); 
uyy = @(x,y,z) -25*u(x,y,z); 
uzz = @(x,y,z) -25*u(x,y,z);
uxy = @(x,y,z) 25*cos(5*x).*cos(5*y).*sin(5*z);
uxz = @(x,y,z) 25*cos(5*x).*sin(5*y).*cos(5*z);
uyz = @(x,y,z) 25*sin(5*x).*cos(5*y).*cos(5*z);

case 6 % sin(8) + sin(12) + sin(14)
u = @(x,y,z) sin(8*x) + sin(12*y) + sin(14*z);
ux = @(x,y,z)  8*cos(8*x);
uy = @(x,y,z) 12*cos(12*y);
uz = @(x,y,z) 14*cos(14*z);
uxx = @(x,y,z) -64*sin(8*x);
uyy = @(x,y,z) -144*sin(12*y);
uzz = @(x,y,z) - 196*sin(14*z);
uxy = @(x,y,z) 0*x.^0; uyx = @(x,y,z) 0*x.^0;
uxz = @(x,y,z) 0*x.^0; uzx = @(x,y,z) 0*x.^0;
uyz = @(x,y,z) 0*x.^0; uzy = @(x,y,z) 0*x.^0;

case 7
a = @(x,y,z) abs(x+y-z-1);
u = @(x,y,z) a(x,y,z).^3;
ux = @(x,y,z) 3*a(x,y,z).^2; uy = @(x,y,z) 3*a(x,y,z).^2; uz = @(x,y,z) 3*a(x,y,z).^2;
uxx = @(x,y,z) 6*a(x,y,z); uxy = @(x,y,z) 6*a(x,y,z); uxz = @(x,y,z) 6*a(x,y,z); 
uyy = @(x,y,z) 6*a(x,y,z); uyz = @(x,y,z) 6*a(x,y,z); uzz = @(x,y,z) 6*a(x,y,z); 

case 8 %homog
u =  @(x,y,z) 40*(x.^2-x).*(y.^2-y).*(z.^2-z);
ux = @(x,y,z)  40*(2*x - 1).*(y.^2-y).*(z.^2-z);
uy = @(x,y,z) 40*(x.^2-x).*(2*y -1).*(z.^2-z);
uz = @(x,y,z) 40*(x.^2-x).*(y.^2-y).*(2*z-1);
uxx = @(x,y,z) 80*(y.^2-y).*(z.^2-z);
uyy = @(x,y,z) 80*(x.^2-x).*(z.^2-z);
uzz = @(x,y,z) 80*(x.^2-x).*(y.^2-y);
uxy = @(x,y,z) 40*(2*x - 1).*(2*y-1).*(z.^2-z);
uxz = @(x,y,z) 40*(2*x - 1).*(y.^2-y).*(2*z-1);
uyz = @(x,y,z) 40*(x.^2-x).*(2*y -1).*(2*z-1);

case 9 % degree 6
u = @(x,y,z) x.* y.^2.*z.^3; 
ux = @(x,y,z) y.^2.*z.^3;
uy = @(x,y,z) 2*x.*y.*z.^3;
uz = @(x,y,z) 3*x.*y.^2.*z.^2;
uxx = @(x,y,z) 0*x.^0;
uyy = @(x,y,z) 2*x.*z.^3;
uzz = @(x,y,z) 6*x.*y.^2.*z;
uxy = @(x,y,z) 2*y.*z.^3;
uxz= @(x,y,z) 3*y.^2.*z.^2;
uyz = @(x,y,z) 6*x.*y.*z.^2;

case 10 % sin(8)*sin(12)*sin(14)
u = @(x,y,z) sin(8*x).*sin(12*y).* sin(14*z);
ux = @(x,y,z) 8*cos(8*x).*sin(12*y).* sin(14*z);
uy = @(x,y,z) 12*sin(8*x).*cos(12*y).* sin(14*z);
uz = @(x,y,z) 14*sin(8*x).*sin(12*y).* cos(14*z);
uxx = @(x,y,z) -64*sin(8*x).*sin(12*y).* sin(14*z);
uyy = @(x,y,z) -144*sin(8*x).*sin(12*y).* sin(14*z);
uzz = @(x,y,z) -196*sin(8*x).*sin(12*y).* sin(14*z);
uxy = @(x,y,z) 96*cos(8*x).*cos(12*y).* sin(14*z);
uyz = @(x,y,z) 168*sin(8*x).*cos(12*y).* cos(14*z);
uxz = @(x,y,z) 112*cos(8*x).*sin(12*y).* cos(14*z);
%a1 = @(x,y,z) 10+cos(x);  a2 = @(x,y,z) exp(y);  a3 = @(x,y,z) 10+sin(x); 

case 11 % sin5x + sin5y + sin5 var coef
u = @(x,y,z) sin(5*x) + sin(5*y) + sin(5*z);
ux = @(x,y,z) 5*cos(5*x); uy = @(x,y,z) 5*cos(5*y); uz = @(x,y,z) 5*cos(5*z);
uxx = @(x,y,z) -25*sin(5*x); uyy = @(x,y,z) -25*sin(5*y); uzz = @(x,y,z) -25*sin(5*z);
uxy = @(x,y,z) 0*x.^0; uyx = @(x,y,z) 0*x.^0;
uxz = @(x,y,z) 0*x.^0; uzx = @(x,y,z) 0*x.^0;
uyz = @(x,y,z) 0*x.^0; uzy = @(x,y,z) 0*x.^0;
a1 = @(x,y,z) 10+cos(x);  a2 = @(x,y,z) exp(y);  a3 = @(x,y,z) 10+sin(x); 
a4 = @(x,y,z) cos(x);  a5 = @(x,y,z) cos(y);  a6 = @(x,y,z) cos(z); 

case 12 % sin5x*sin5y*sin5z -- coef (1,1,1,1,1,1)
u= @(x,y,z) sin(5*x).*sin(5*y).*sin(5*z);
ux = @(x,y,z) 5*cos(5*x).*sin(5*y).*sin(5*z);
uy = @(x,y,z) 5*sin(5*x).*cos(5*y).*sin(5*z);
uz = @(x,y,z) 5*sin(5*x).*sin(5*y).*cos(5*z);
uxx = @(x,y,z) -25*u(x,y,z);
uyy = @(x,y,z) -25*u(x,y,z);
uzz = @(x,y,z) -25*u(x,y,z);
uxy = @(x,y,z) 25*cos(5*x).*cos(5*y).*sin(5*z);
uxz = @(x,y,z) 25*cos(5*x).*sin(5*y).*cos(5*z);
uyz = @(x,y,z) 25*sin(5*x).*cos(5*y).*cos(5*z);
a1 = @(x,y,z) x.^0;  a2 = @(x,y,z) y.^0;  a3 = @(x,y,z) z.^0; 
a4 = @(x,y,z) x.^0;  a5 = @(x,y,z) y.^0;  a6 = @(x,y,z) z.^0; 

case 13 % sin5x*sin5y*sin5z -- coef (1,1,1,10,10,10)
u= @(x,y,z) sin(5*x).*sin(5*y).*sin(5*z);
ux = @(x,y,z) 5*cos(5*x).*sin(5*y).*sin(5*z);
uy = @(x,y,z) 5*sin(5*x).*cos(5*y).*sin(5*z);
uz = @(x,y,z) 5*sin(5*x).*sin(5*y).*cos(5*z);
uxx = @(x,y,z) -25*u(x,y,z); 
uyy = @(x,y,z) -25*u(x,y,z); 
uzz = @(x,y,z) -25*u(x,y,z);
uxy = @(x,y,z) 25*cos(5*x).*cos(5*y).*sin(5*z);
uxz = @(x,y,z) 25*cos(5*x).*sin(5*y).*cos(5*z);
uyz = @(x,y,z) 25*sin(5*x).*cos(5*y).*cos(5*z);
a1 = @(x,y,z) x.^0;  a2 = @(x,y,z) y.^0;  a3 = @(x,y,z) z.^0; 
a4 = @(x,y,z) 10*x.^0;  a5 = @(x,y,z) 10*y.^0;  a6 = @(x,y,z) 10*z.^0; 

case 14 % sin(8) + sin(12) + sin(14)
u = @(x,y,z) sin(8*x) + sin(12*y) + sin(14*z);
ux = @(x,y,z)  8*cos(8*x);
uy = @(x,y,z) 12*cos(12*y);
uz = @(x,y,z) 14*cos(14*z);
uxx = @(x,y,z) -64*sin(8*x);
uyy = @(x,y,z) -144*sin(12*y);
uzz = @(x,y,z) - 196*sin(14*z);
uxy = @(x,y,z) 0*x.^0; uyx = @(x,y,z) 0*x.^0;
uxz = @(x,y,z) 0*x.^0; uzx = @(x,y,z) 0*x.^0;
uyz = @(x,y,z) 0*x.^0; uzy = @(x,y,z) 0*x.^0;
a1 = @(x,y,z) 10+cos(x);  a2 = @(x,y,z) exp(y);  a3 = @(x,y,z) 10+sin(x);


end;

g = @(x,y,z) u(x,y,z) ;
f = @(x,y,z) a1(x,y,z).*uxx(x,y,z) + a2(x,y,z).*uyy(x,y,z)  + ...
   a3(x,y,z).*uzz(x,y,z) + 2*a4(x,y,z).*uxy(x,y,z) + 2*a5(x,y,z).*uxz(x,y,z) + ...
      2*a6(x,y,z).*uyz(x,y,z);


