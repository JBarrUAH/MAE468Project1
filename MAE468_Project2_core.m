% This is the core code file for the MAE468 Project 2 submission
% The team consists of Joseph Barragree, Sarah Polickoski, Micajah
% Schweikert, and Stephen Ward.

%% Notes
%source used for some constants: http://www.dept.aoe.vt.edu/~lutze/AOE2104/consts.pdf
%source for Mars and Earth patch-conic https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html

%% Housekeeping
% Run to remove figures, workspace variables and command window content
format compact
close all
clear
clc

%% Orbital Variable Initialization
% Includes orbital elements, locations, constants, anonymous functions, etc

mu=1; %canonical mu, AU^3/TU^2 heliocentric or DU^3/TU^2 geocentric
% a,e,i,OMEGA,omega,theta (AU,unitless,deg,deg,deg,deg)
oeE=[1.000000,0.01671,0.00005,-11.26064,114.20783,-2.48284]; %Earth
oeM=[1.523662,0.093412,1.85061,49.57854,286.4623,19.41248]; %Mars
% a,e,i,O,w,th is naming convention used in functions
t0=datetime(2000,1,1,11,58,0); %setting initial time to the J2000 parameter
ToFfun=@(tt) etime(datevec(tt),datevec(t0))/5.0226757e6; %anonymous function for finding ToF in AU

zg=18; %initial guess for z (for Gauss Orbit). Should be within the range of +-(2pi)^2 
% NOTE: Use 18 if elliptical, 0 if parabolic, -18 if hyperbolic

vCAtoM=@(v) v*149597870.7/5.0226757e6; %km/s (AU/TU sun to km/s)

%% Orbital Elements to Initial Vectors
% Finds perifocal vectors from orbital elements, then converts to
% vectors in heliocentric frame. Setup work for other tasks.
[rxyzE0,vxyzE0]=OEtoXYZ(oeE(1),oeE(2),oeE(3),oeE(4),oeE(5),oeE(6),1);
[rxyzM0,vxyzM0]=OEtoXYZ(oeM(1),oeM(2),oeM(3),oeM(4),oeM(5),oeM(6),1);
%% Tasks a and b
% Obtains positions and velocities of Earth and Mars with a 190 day
% interval launching on 17 Sep 2022
fprintf("---TASK a---\n");
tE1=datetime(2022,9,17,15,0,0); %setting departure date, 17 Sep 2022 at 1500hrs UTC was determined to be optimal through manual iteration
tM2=tE1+days(190); %adding 190 days to find the future Mars position
disp("Earth launch date: "+char(tE1));
disp("Mars arrival date: "+char(tM2));

fprintf("\n---TASK b---\n");
[rxyzE1,vxyzE1]=uToF(rxyzE0,vxyzE0,ToFfun(tE1),1); %obtaining Earth vectors at departure
[rxyzM2,vxyzM2]=uToF(rxyzM0,vxyzM0,ToFfun(tM2),1); %obtaining Mars vectors at arrival
fprintf("Earth at Departure\n\t Position: %5.4f, %5.4f, %5.4f AU\n\t Velocity: %5.4f, %5.4f, %5.4f AU/TU",rxyzE1(1),rxyzE1(2),rxyzE1(3),vxyzE1(1),vxyzE1(2),vxyzE1(3));
fprintf("\nMars at Arrival\n\t Position: %5.4f, %5.4f, %5.4f AU\n\t Velocity: %5.4f, %5.4f, %5.4f AU/TU\n",rxyzM2(1),rxyzM2(2),rxyzM2(3),vxyzM2(1),vxyzM2(2),vxyzM2(3));

%% Task c
% Spacecraft interplanetary transfer parameters
fprintf("\n---TASK c---\n");
[vxyz1,vxyz2]=Gorb(rxyzE1,rxyzM2,ToFfun(t0+days(190)),1,zg); %spacecraft heliocentric departure and arrival velocities
[oeSC1(1),oeSC1(2),oeSC1(3),oeSC1(4),oeSC1(5),oeSC1(6)]=XYZtoOE(rxyzE1,vxyz1,1); %Spacecraft departure orbital elements
[~,~,~,~,~,oeSC2]=XYZtoOE(rxyzM2,vxyz2,1); %Spacecraft arrival true anomaly
fprintf("Spacecraft RoI Orbital Elements\n\t Eccentricity: %5.3f\n\t Semi-major axis: %5.3f AU\n\t Inclination: %5.3f degrees\n\t RAAN: %5.3f degrees\n\t Argument of periapsis: %5.3f degrees\n\t Departure true anomaly: %5.3f degrees\n\t Arrival true anomaly: %5.3f degrees\n",oeSC1(1),oeSC1(2),oeSC1(3),oeSC1(4),oeSC1(5),oeSC1(6),oeSC2);
fprintf("Spacecraft at Earth RoI departure\n\t Heliocentric Velocity: %5.4f, %5.4f, %5.4f AU/TU",vxyz1(1),vxyz1(2),vxyz1(3));
fprintf("\nSpacecraft at Mars RoI Arrival\n\t Heliocentric Velocity: %5.4f, %5.4f, %5.4f AU/TU\n",vxyz2(1),vxyz2(2),vxyz2(3));

%% Task d
% find delta v and calculations for the trajectory assuming ???
% launch and 400km circular Mars parking orbit
fprintf("\n---TASK d---");
vinf1=vCAtoM(norm(vxyz1-vxyzE1)); %finding vinf at Earth in km/s
vinf2=vCAtoM(norm(vxyz2-vxyzM2)); %finding vinf at Mars in km/s
dvE=sqrt(2*(vinf1^2/2+3.986e5/6578.1))-sqrt(3.986e5/6578.1); %finding dv for Earth orbit
dvM=sqrt(4.2828e4/3796.2)-sqrt(2*(vinf2^2/2+4.2828e4/3796.2));%finding dv for Mars 400km altitude parking orbit
fprintf("\n-OUTDATED- Earth parking orbit to vinf dv: %5.4f km/s\n",abs(dvE)); %displaying delta v
fprintf("Mars vinf to parking orbit dv: %5.4f km/s\n",abs(dvM)); %displaying delta v to get Mars parking orbit

%% Task e
% Size systems and determine costs
fprintf("\n---TASK e---\n");
CasIns=[11.83,14.46,32;55.9,57.8,365;3.1,3,3.6;27.7,9.25,1.5]; %Cassini Instrument parameters [W,kg,max kb/s]



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
rpqw=[rm*cosd(th);rm*sind(th);0]; %perifocal vectors
vpqw=sqrt(mu/p)*[-sind(th);e+cosd(th);0];
R=[cosd(O)*cosd(w)-sind(O)*sind(w)*cosd(i),-cosd(O)*sind(w)-sind(O)*cosd(w)*cosd(i),sind(O)*sind(i);sind(O)*cosd(w)+cosd(O)*sind(w)*cosd(i),-sind(O)*sind(w)+cosd(O)*cosd(w)*cosd(i),-cosd(O)*cosd(i);sind(w)*sind(i),cosd(w)*sind(i),cosd(i)];
rxyz=R*rpqw; %heliocentric-ecliptic vectors
vxyz=R*vpqw;
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

function [e,a,i,O,w,th] = XYZtoOE(rxyz,vxyz,mu)
% Heliocentric vectors to orbital elements
% inputs rxyz, vxyz, mu
hxyz=cross(rxyz,vxyz);
hm=norm(hxyz);
n=cross([0;0;1],hxyz);
nm=norm(n);
rm=norm(rxyz);
ev=1/mu*((dot(vxyz,vxyz)-mu/rm)*rxyz-(dot(rxyz,vxyz)*vxyz));

e=norm(ev);
p=hm^2/mu;
a=p/(1-e^2);
i=acosd(dot([0;0;1],hxyz)/hm);
O=acosd(n(1)/nm);
if n(2)<0
    O=360-O;
end
w=acosd(dot(n,ev)/(nm*e));
if ev(3)<0
    w=360-w;
end
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

