% This is the core code file for the MAE468 Project 1 submission
% The team consists of Joseph Barragree, Sarah Polickoski, Micajah
% Schweikert, and Stephen Ward.

%% Notes
% This section is for leaving note comments

%source used for some constants: http://www.dept.aoe.vt.edu/~lutze/AOE2104/consts.pdf
%source for Mars and Earth patch-conic https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html

%% Housekeeping
% Run to remove figures, workspace variables and command window content
format compact
close all
clear
clc

%% Variable Initialization
% Initial setup. Includes orbital elements, locations, constants, etc

mu=1; %canonical mu, AU^3/TU^2 heliocentric or DU^3/TU^2 geocentric
% a,e,i,OMEGA,omega,theta (AU,unitless,deg,deg,deg,deg)
oeE=[1.000000,0.01671,0.00005,-11.26064,114.20783,-2.48284]; %Earth
oeM=[1.523662,0.093412,1.85061,49.57854,286.4623,19.41248]; %Mars
oeJ=[5.203363,0.048393,1.3053,100.55615,-85.8023,19.55053]; %Jupiter
% a,e,i,O,w,th is naming convention used in functions
t0=datetime(2000,1,1,11,58,0); %setting initial time to the J2000 parameter
ToFfun=@(tt) etime(datevec(tt),datevec(t0))/5.0226757e6; %anonymous function for finding ToF in AU, might use actual function

zg=18; %initial guess for z (for Gauss Orbit). Should be within the range of +-(2pi)^2 
% NOTE: Use 18 if elliptical, 0 if parabolic, -18 if hyperbolic

vCAtoM=@(v) v*149597870.7/5022604.8; %km/s (AU/TU sun to km/s)

%% Orbital Elements to Initial Vectors
% Finds perifocal vectors from orbital elements, then converts to
% vectors in heliocentric frame. Setup work for other tasks.
[rxyzE0,vxyzE0]=OEtoXYZ(oeE(1),oeE(2),oeE(3),oeE(4),oeE(5),oeE(6),1);
[rxyzM0,vxyzM0]=OEtoXYZ(oeM(1),oeM(2),oeM(3),oeM(4),oeM(5),oeM(6),1);
[rxyzJ0,vxyzJ0]=OEtoXYZ(oeJ(1),oeJ(2),oeJ(3),oeJ(4),oeJ(5),oeJ(6),1);

%% Task 1
% Obtains rxyz and vxyz vectors for the planets on Dec 25, 2025 at 0837 UTC
% as well as the true anomalies
t1=datetime(2025,12,25,08,37,00); %specified time as datetime vector
[rxyzE01,vxyzE01]=uToF(rxyzE0,vxyzE0,ToFfun(t1),1); %obtaining final vectors for Earth
[thE01]=Tanomaly(rxyzE01,vxyzE01,1); %Earth true anomaly
[rxyzM01,vxyzM01]=uToF(rxyzM0,vxyzM0,ToFfun(t1),1); %obtaining final vectors for Mars
[thM01]=Tanomaly(rxyzM01,vxyzM01,1);
[rxyzJ01,vxyzJ01]=uToF(rxyzJ0,vxyzJ0,ToFfun(t1),1); %obtaining final vectors for Jupiter
[thJ01]=Tanomaly(rxyzJ01,vxyzJ01,1);
fprintf("Earth Data\n\t Position: %5.4f %5.4f %5.4f AU\n\t Velocity: %5.4f %5.4f %5.4f AU/TU\n\t True anomaly: %5.2f degrees",rxyzE01(1),rxyzE01(2),rxyzE01(3),vxyzE01(1),vxyzE01(2),vxyzE01(3),thE01); %displaying results
fprintf("\nMars Data\n\t Position: %5.4f %5.4f %5.4f AU\n\t Velocity: %5.4f %5.4f %5.4f AU/TU\n\t True anomaly: %5.2f degrees",rxyzM01(1),rxyzM01(2),rxyzM01(3),vxyzM01(1),vxyzM01(2),vxyzM01(3),thM01); %displaying results
fprintf("\nJupiter Data\n\t Position: %5.4f %5.4f %5.4f AU\n\t Velocity: %5.4f %5.4f %5.4f AU/TU\n\t True anomaly: %5.2f degrees\n\n",rxyzJ01(1),rxyzJ01(2),rxyzJ01(3),vxyzJ01(1),vxyzJ01(2),vxyzJ01(3),thJ01); %displaying results

%% Task 2
% Obtains positions of Earth and Mars with a 190 day interval, then
% calculates the required V0 of the spacecraft from Earth's R0 to make
% the 190 day shortway transfer. Arbitrarily set the Earth date, needs to be between
% 2021 and 2030.

tE1=datetime(2022,9,12,13,0,0); %setting departure date, 12 Sep 2022 at 1300hrs UTC was determined to be optimal through manual iteration
tM2=tE1+days(190); %adding 190 days to find the future Mars position

% Departure
[rxyzE1,vxyzE1]=uToF(rxyzE0,vxyzE0,ToFfun(tE1),1); %obtaining Earth vectors at departure
[thE1]=Tanomaly(rxyzE1,vxyzE1,1); %Earth true anomaly
[rxyzM1,vxyzM1]=uToF(rxyzM0,vxyzM0,ToFfun(tE1),1); %obtaining Mars vectors at arrival
[thM1]=Tanomaly(rxyzM1,vxyzM1,1); %Mars true anomaly
% Arrival
[rxyzE2,vxyzE2]=uToF(rxyzE0,vxyzE0,ToFfun(tM2),1); %obtaining Earth vectors at departure
[thE2]=Tanomaly(rxyzE2,vxyzE2,1); %Earth true anomaly
[rxyzM2,vxyzM2]=uToF(rxyzM0,vxyzM0,ToFfun(tM2),1); %obtaining Mars vectors at arrival
[thM2]=Tanomaly(rxyzM2,vxyzM2,1); %Mars true anomaly
disp("Soonest ideal launch date: "+char(tE1));
fprintf("Earth at Departure\n\t Position: %5.4f %5.4f %5.4f AU\n\t Velocity: %5.4f %5.4f %5.4f AU/TU\n\t True anomaly: %5.2f degrees",rxyzE1(1),rxyzE1(2),rxyzE1(3),vxyzE1(1),vxyzE1(2),vxyzE1(3),thE1); %displaying results
fprintf("\nMars at Arrival\n\t Position: %5.4f %5.4f %5.4f AU\n\t Velocity: %5.4f %5.4f %5.4f AU/TU\n\t True anomaly: %5.2f degrees\n",rxyzM2(1),rxyzM2(2),rxyzM2(3),vxyzM2(1),vxyzM2(2),vxyzM2(3),thM2);

S = [0;0];
E1 = rxyzE1(1:2);
M1 = rxyzM1(1:2);
E2 = rxyzE2(1:2);
M2 = rxyzM2(1:2);
Rfactor = [1 .5 .25];
Rsun = 0.2;
lbl = ["Sun";"Earth";"Mars"];
clr = ["y";"g";"r"];

%% Task 3
% Determine net deltaV assuming departure from 200km Earth circ parking orbit
% and arrival at 1000 km Mars circ park orb
[vxyz1,vxyz2]=Gorb(rxyzE1,rxyzM2,ToFfun(t0+days(190)),1,zg); %spacecraft heliocentric departure and arrival velocities
vinf1=vCAtoM(norm(vxyz1-vxyzE1)); %finding vinf at Earth in km/s
vinf2=vCAtoM(norm(vxyz2-vxyzM2)); %finding vinf at Mars in km/s
dvE=sqrt(2*(vinf1^2/2+3.986e5/6578.1))-sqrt(3.986e5/6578.1); %finding dv at Earth orbit
dvM=sqrt(4.2828e4/4396.2)-sqrt(2*(vinf2^2/2+4.2828e4/4396.2));%finding dv at Mars orbit
fprintf("\nNet mission dv: %5.4f km/s\n",abs(dvE)+abs(dvM)); %displaying net delta v for entire mission

%% Trajectory plotting
[e,a,Omega,omega,Inc,theta,h] = OrbitalElements(rxyzE1,vxyz1,1);
figure
hold on
plotBody(Rsun*Rfactor(1),S(1),S(2),lbl(1),clr(1),.1)
plotBody(Rsun*Rfactor(2),E1(1),E1(2),lbl(2),clr(2),.1)
plotBody(Rsun*Rfactor(3),M1(1),M1(2),lbl(3),clr(3),.1)

plotBody(Rsun*Rfactor(1),S(1),S(2),lbl(1),clr(1))
plotBody(Rsun*Rfactor(2),E2(1),E2(2),lbl(2),clr(2))
plotBody(Rsun*Rfactor(3),M2(1),M2(2),lbl(3),clr(3))
% Not sure if the theta is correct here
plotOrbit(a,e,-theta,'b.')

title({"Trajectory Patch Conic";"Lower opacity is initial location and higher opacity is final location"})
xlabel("x Space dimension [Au]")
ylabel("y Space dimension [Au]")
hold off

%% Functions
% Organized here for ease of editing
function [rxyz,vxyz] = OEtoXYZ(a,e,i,O,w,th,mu)
% Convertes orbital parameters to heliocentric vectors
% Inputs
%   a - semimajor axis [Au]
%   e - eccentricity []
%   i - inclination angle [deg]
%   O - Longitude of Ascending node [deg]
%   w - Argument of Periapsis [deg]
%   th - true anomaly [deg]
%   mu - gravitional parameter [Au^3/Tu^2] assumed 1
% Outputs
%   rxyz - radius vector of orbiter [Au]
%   vxyz - velocity vector of orbiter [Au/Tu]

p=a*(1-e^2);
rm=(p)/(1+e*cosd(th));
rpkw=[rm*cosd(th);rm*sind(th);0]; %perifocal vectors
vpkw=sqrt(mu/p)*[-sind(th);e+cosd(th);0];
R=[cosd(O)*cosd(w)-sind(O)*sind(w)*cosd(i),-cosd(O)*sind(w)-sind(O)*cosd(w)*cosd(i),sind(O)*sind(i);sind(O)*cosd(w)+cosd(O)*sind(w)*cosd(i),-sind(O)*sind(w)+cosd(O)*cosd(w)*cosd(i),-cosd(O)*cosd(i);sind(w)*sind(i),cosd(w)*sind(i),cosd(i)];
rxyz=R*rpkw; %heliocentric-ecliptic vectors
vxyz=R*vpkw;
end

function [r1,v1] = uToF(r0,v0,ToF,mu)
% Universal Time of Flight calculator propogates r0 and v0 through TOF to r1 v1
% Inputs
%   r0 - initial radius vector [Au]
%   v0 - initial velocity vector [Au/Tu]
%   TOF - time of flight in Tu
%   mu - gravitational parameter [Au^3/Tu^2] assumed 1
% Outputs
%   r1 - final radius vector [Au]
%   v1 - final velocity vector [Au/Tu]
% Universal ToF using z (for all orbit types)
r0m=norm(r0); %radius magnitude
v0m=norm(v0); %velocity magnitude
ene=v0m^2/2-mu/r0m; %orbital energy
a=-mu/(2*ene); % semi-major axis
xo=ToF*sqrt(mu)/a; %initial guess for x
zo=xo^2/a; %initial guess for z
%Selecting anonymous functions used during iteration based on orbit
if ene == 0 %parabolic parameters
    S=@(z) 1/factorial(3)-z/factorial(5)+z^2/factorial(7)-z^3/factorial(9);
    C=@(z) 1/factorial(2)-z/factorial(4)+z^2/factorial(6)-z^3/factorial(8);
elseif zo > 0 %elliptical parameters
    S=@(z) (sqrt(z)-sin(sqrt(z)))/sqrt(z^3);
    C=@(z) (1-cos(sqrt(z)))/z;
elseif zo < 0 %hyperbolic parameters
    S=@(z) (sinh(sqrt(-z))-sqrt(-z))/sqrt((-z)^3);
    C=@(z) (1-cosh(sqrt(-z)))/z;
else
    fprintf("The universe has encountered a fatal error. Your spacecraft will now self-destruct");
end
%iterative solver
while 1
    t=((xo^3*S(zo))+(dot(r0,v0)/sqrt(mu)*xo^2*C(zo))+(r0m*xo*(1-zo*S(zo))))/sqrt(mu);
    rn=xo^2*C(zo)+dot(r0,v0)/sqrt(mu)*xo*(1-zo*S(zo))+r0m*(1-zo*C(zo));
    if abs((ToF-t)/ToF) <= 1e-5 %terminating if difference between times is small enough
        x=xo;
        break
    end
    xo=xo+(ToF-t)/(rn/sqrt(mu)); %newton-step to next x guess
    zo=xo^2/a;
end
%processing the solution
f=1-a/r0m*(1-cos(x/sqrt(a)));
g=a^2/sqrt(mu*a)*((dot(r0,v0)/sqrt(mu*a)*(1-cos(x/sqrt(a))))+(r0m/a*sin(x/sqrt(a))));
r1=f*r0+g*v0; %final position vector
r1m=norm(r1);
df=-sqrt(mu*a)/(r0m*r1m)*sin(x/sqrt(a));
dg=1-a/r1m+a/r1m*cos(x/sqrt(a));
v1=df*r0+dg*v0; %final velocity vector
end

function [th] = Tanomaly(rxyz,vxyz,mu)
% Calculates the true anomaly of the orbit with these vectors
% Inputs
%   rxyz - radius vector [Au]
%   vxyz - velocity vector [Au/Tu]
%   mu - gravitational parameter [Au^2/TU^2] assumed 1
% Outputs
%   th - True anomaly location [deg]
rm=norm(rxyz);
ev=1/mu*((dot(vxyz,vxyz)-mu/rm)*rxyz-(dot(rxyz,vxyz)*vxyz));
e=norm(ev);
th=acosd(dot(ev,rxyz)/(e*rm));
if dot(rxyz,vxyz)<0
    th=360-th;
end
end

function [v0s,v1s] = Gorb(r0,r1,ToF,mu,zg)
% Gauss Orbit solver (all types of orbit)
% Inputs 
%   r0 - initial radius vector [Au]
%   ToF - Time of flight [TU] 
%   mu - gravitational constant [Au^2/Tu^2] assumed 1
%   r1 - final radius vector [Au] 
%   zg - initial Guess Use 18 if elliptical, 0 if parabolic, -18 if
%   hyperbolic
% Outputs
%   v0s - shortpath initial velocity vector
%   v1s - shortpath final velocity vector
r0m=norm(r0); %initial radius magnitude
r1m=norm(r1); %final radius magnitude
Dths=acosd(dot(r0,r1)/(r0m*r1m)); %delta theta (short path), in degrees

As=sqrt(r0m*r1m)*sind(Dths)/sqrt(1-cosd(Dths)); %short path A
%Selecting anonymous functions used during iteration based on orbit
if zg == 0 %parabolic parameters
    S=@(z) 1/factorial(3)-z/factorial(5)+z^2/factorial(7)-z^3/factorial(9);
    C=@(z) 1/factorial(2)-z/factorial(4)+z^2/factorial(6)-z^3/factorial(8);
elseif zg > 0 %elliptical parameters
    S=@(z) (sqrt(z)-sin(sqrt(z)))/sqrt(z^3);
    C=@(z) (1-cos(sqrt(z)))/z;
elseif zg < 0 %hyperbolic parameters
    S=@(z) (sinh(sqrt(-z))-sqrt(-z))/sqrt((-z)^3);
    C=@(z) (1-cosh(sqrt(-z)))/z;
else
    fprintf("Invalid z guess");
end
%iterating to find parameters
while 1
    yo=r0m+r1m-As*(1-zg*S(zg))/sqrt(C(zg));
    xo=sqrt(yo/C(zg));
    ts=(xo^3*S(zg)+As*sqrt(yo))/sqrt(mu);
    if abs((ToF-ts)/ToF) < 1e-5 %convergence when relative error is small enough
        xs=xo;%outputting short path parameters for processing
        zs=zg;
        break
    end
    zg=zg+(ToF-ts)/((xo^3*((C(zg)-3*S(zg))/(2*zg)-(3*S(zg)*((1-zg*S(zg)-2*C(zg))/(2*zg)))/(2*C(zg)))+As/8*(3*S(zg)*sqrt(yo)/C(zg)+As/xo))/sqrt(mu));
end
%processing results
fs=1-xs^2/r0m*C(zs);
gs=ts-xs^3/sqrt(mu)*S(zs);
dgs=1-xs^2/r1m*C(zs);
v0s=(r1-fs*r0)/gs; %outputting short path initial and final velocity vectors
v1s=(dgs*r1-r0)/gs;
end

function plotBody(R,xdiff,ydiff,lbl,clr,alf,SOI)
% Plots a celestial body at location xdiff ydiff with pseudo radius R [Au]
% with label  and optional Sphere of influence
% Input
%   R - psuedo radius (aparant radius) [Au]
%   xdiff - shifted x location [Au]
%   ydiff - shifted y location [Au]
%   clr - color of body [string]
%   SOI - sphere of influence
%   alf - alpha of plot (float)
% Output
%   plot

SunR = R;
xsun = linspace(-SunR,SunR,100);
ysunp = sqrt(SunR^2-xsun.^2);
ysunn = -ysunp;
xsun = [xsun fliplr(xsun)]+xdiff;
ysun = [ysunp fliplr(ysunn)]+ydiff;

if or(xdiff~=0,ydiff~=0)
    R = sqrt(xdiff^2+ydiff^2);
    x = linspace(-R,R,100);
    yp = sqrt(R^2-x.^2);
    yn = -sqrt(R^2-x.^2);
    x = [x fliplr(x)];
    y = [yp fliplr(yn)];
else
    x = 0;
    y = 0;
end

if nargin==7 
    t = linspace(0,2*pi,100);
    r = ones(size(t))*SOI;
    [xSOI,ySOI] = pol2cart(t,r);
     xSOI = xSOI+xdiff;
     ySOI = ySOI+ydiff;
    plot(xSOI,ySOI,'b.')
else
    xSOI=xdiff;
    ySOI = ydiff;
end
a = plot(xsun,ysun,clr,x,y,'k',xSOI,ySOI,'b.');
if exist('alf')
    fill(xsun,ysun,clr,'FaceAlpha',alf);
    a(1).Color(4) = alf;
else
    fill(xsun,ysun,clr,'FaceAlpha',1);
    a(1).Color(4) = 1;
end
text(mean(xsun),mean(ysun),lbl)
axis('equal')
end

function [e,a,Omega,omega,Inc,theta,h] = OrbitalElements(r,v,mue)
% Calculates the Orbital elements of a given orbit with r,v
% Requires:
%   r - distance column vector of the spacecraft
%   v - velocity column vector of the spacecraft
% Returns:
%   e - ecentricity 
%   a - semi-major axis
%   Omega - Longitude of ascending node [degrees];
%   omega - argument of periapsis [degrees];
%   Inc - inclination [degrees]
%   theta - true anaomaly [degrees]
%   h - angular momentum

K = [0;0;1];
J = [0;1;0];
I = [1;0;0];

h = cross(r,v);
n = cross(K,h);
e =1/mue*((norm(v)^2-mue/norm(r))*r-(dot(r,v))*v);

a = (norm(h)^2/mue)/(1-norm(e)^2);

i = acosd(dot(K,h)/norm(h));

if dot(n,J)>0
    Omega = acosd(dot(I,n)/norm(n));
else
    Omega = 360 - acosd(dot(I,n)/norm(n));
end

if dot(e,K)>0
    omega = acosd(dot(n,e)/(norm(n)*norm(e)));
else
    omega = 360 - acosd(dot(n,e)/(norm(n)*norm(e)));
end

if dot(r,v)>0
    theta = acosd(dot(e,r)/(norm(e)*norm(r)));
else
    theta = 360 -acosd(dot(e,r)/(norm(e)*norm(r)));
end

e = norm(e);
h = norm(h);
Inc = i;
% info about orbit
if i>90
    grad = "retrograde";
elseif i==90
    grad = "polar";
else
    grad = "prograde";
end

if e==0
    typ = "circular";
elseif e<1;
    typ = "elliptical";
elseif e==1;
    typ = "parabolic";
else
    typ = "hyperbolic";

if isnan(Omega)
    k = "equatorial";
else
    k = "";
end

fprintf(" Orbit is "+k+" "+grad+" and is "+typ+"\n")
end
end

function plotOrbit(a,e,theta,clr)
    if nargin<4
        clr = 'k';
    end
    if nargin<3
        theta = 0;
    end
    
    R = [cosd(theta) -sind(theta);
        sind(theta) cosd(theta)];
    t = linspace(0,2*pi,100);
    p = a.*(1-e.^2);
    r = p./(1+e*cos(t));
    [x1,y1] = pol2cart(t,r);
    x =[];
    y = [];
    for i = 1:length(x1)
        X = R*[x1(i);y1(i)];
        x(i) = X(1);
        y(i) = X(2);
    end
    plot(x,y,clr)
    axis('equal')
end