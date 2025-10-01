% T3flam -- test lambda for IPBF
% 2025: 7/30
lambda= input('input lambda ');  lams = lambda;

for i = 1:4
  [c,G,t1,t2] = ipbf0d3(d,x,y,z,TET,E,F,evertex,fvertex,...
    a1,a2,a3,a4,a5,a6,f,g,xb,yb,zb,S,lambda,lams);
 val = valsp3(d,x,y,z,TET,E,F,evertex,fvertex,xr,yr,zr,c);
 err = u(xr,yr,zr) - val; em(i) = norm(err,inf); er(i) = erms(err);
  fprintf(' %.3f & %.2f & %.2f & %d & %d & %.2e & %.2e & %.2e \\cr \n',...
      lambda,t1,t2,length(TET),length(c),condest(G), norm(err,inf),erms(err));
 lambda = 10*lambda; lams = 10*lams;
end;
