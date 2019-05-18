function Dy = Reentrypr(t,y)
global mu RE CD CL S m B E 
%global k0 omega0

r = y(1);
la = y(2);
L = y(3);
psi = y(4);
ga = y(5);
V = y(6);


h = r-RE;
rho0 = 1.22e+09;            % air density kg/Km^3
beta = 0.1414;              % Scale Factor 1/km
rho = rho0*exp(-beta*h);
B = 0.5*CD*S/m;             % Ballistic Factor
E = CL/CD;                  % Efficency

%mi_star=k0*sin(omega0*t);
mi_star=0;
%mi_star=-pi/2;
%mi_star=pi/2;

a1 = V*sin(ga);
a2 = (V*cos(ga)*cos(psi))/(r*cos(L));
a3 = -(V*cos(ga)*sin(psi))/r;

a4 = -((V*cos(ga)*cos(psi)*tan(L))/r)+rho*B*E*V*sin(mi_star);
a5 = -mu/(V*r^2)*cos(ga)+V/r*cos(ga)+rho*B*E*V*cos(mi_star);
a6 = -mu/r^2*sin(ga)- rho*B*V^2; 


Dy =[a1 a2 a3 a4 a5 a6]';

end

