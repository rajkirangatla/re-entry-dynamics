clear all
close all
clc

global mu RE CD CL S m B E k0 omega0

mu = 398601;            % km^3/s^2
RE = 6378;              % Km
CD = 1.34;              % Soyuz Data: CD = 1.341 
CL = 0.35;              % Soyuz Data: CL = 0.349
S = 3.8;                % m^2  (Soyuz Data)
S = S*1e-06;
m = 2400;               % Kg (Soyuz Data)
mi_star = 0;
k0=0.;
Tp=20;
omega0=2*pi/Tp;

%  initial conditions

V0 = 7.9;               % Km/s
gamma0 = -5.5;         % deg
lambda0 = 45;           % deg
L0 = 45;                % deg
psi0 = 28;              % deg
h0 = 300;               % Km
rho0 = 1.22e+09;        % air density kg/Km^3
beta = 0.1414;          % Scale Factor 1/km

r0 = h0+RE;
B = 0.5*CD*S/m;         % Ballistic Factor
E = CL/CD;              % Efficency

Y0=[r0 lambda0*pi/180 L0*pi/180 psi0*pi/180 gamma0*pi/180 V0 ]  % initial values for the reentry vehicles
tspan = [0 10000];
Tol =1e-12;
Tol0 = 1e-12;
options=odeset('RelTol',Tol0,'AbsTol',Tol,'events',@events3);
[t,y] = ode45(@Reentrypr, tspan, Y0, options);

xz = length(t);

R=y(:,1);
h = R - RE;
rho = rho0*exp(-beta*h);
La=y(:,2);
L=y(:,3);
psi=y(:,4);
gama=y(:,5);
Vz=y(:,6);


X = R.*cos(L).*cos(La);
Y = R.*cos(L).*sin(La);
Z = R.*sin(L);

Lf = L(xz);                     % Final Latitude 36.87 in degrees
laf = La(xz);                   % Final longitude (lambda)  89.145 in degrees 

[xe ye ze] = sphere;
XE = RE*xe;
YE = RE*ye;
ZE = RE*ze;

i = acos(cos(psi).*cos(L));  % inclination cos(i)=cos(psi)cos(L)
theta = asin(sin(L)./sin(i)); % theta sin(i)=sin(theta)*sin(i)
Ohm = asin((sin(psi).*sin(La)-cos(psi).*sin(L).*cos(La))./sin(i)); % omega
Energy = (Vz.^2/2)-(mu./R);
a = -mu/2*Energy;                               % semiaxis
e = sqrt(1-(R.*Vz.*cos(gama)).^2./(mu.*a));        % eccentricity


%----------------------------- Figures --------------------------------%

figure(1)
subplot(2,1,1)
plot(t,Vz)
title('Reentry 3D')
ylabel('Velocity (Km/s)')
xlabel('time (s)')

subplot(2,1,2)
plot(h,Vz)
title('Reentry 3D')
xlabel('Altitude (km)')
ylabel('Velocity (Km/s)')

figure(2)
plot(t,h)
title('Reentry 3D')
ylabel('Altitude (Km)')
xlabel('time (s)')

figure(3)
plot(t,gama*180/pi)
title('Reentry 3D')
ylabel('Flight Path Angle (deg)')
xlabel('time (s)')

figure(4)
surf(XE,YE,ZE)
hold on
earth_g=imread('world3.jpg');
imshow(earth_g)
hold on
warp(-XE,YE,-ZE,earth_g)
plot3(X,Y,Z)
% hold on
% surf(XE,YE,ZE)
title('Trajectory')

disp('Final Latitude and Longitude (deg)')
Latitude = Lf*180/pi
Longitude = laf*180/pi

   