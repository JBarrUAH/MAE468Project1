% This is the core code file for the MAE468 Project 1 submission
% The team consists of Joseph Barragree, Sarah Polickoski, Micajah
% Schweikert, and Stephen Ward.

% This code requires ___ to run. (put.m file names here if others are needed)
%% Notes
% This section is for leaving note comments

%source used for some constants: http://www.dept.aoe.vt.edu/~lutze/AOE2104/consts.pdf

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
t0=[2000,1,1,11,58,0]; %setting initial time to the J2000 parameter
ToFfun=@(tt) etime(tt,t0)/5.0226757e6; %anonymous function for finding ToF in AU, might use actual function

zg=18; %initial guess for z (for Gauss Orbit). Should be within the range of +-(2pi)^2 
% NOTE: Use 18 if elliptical, 0 if parabolic, -18 if hyperbolic

%% Orbital Elements to Initial Vectors
% Finds perifocal vectors from orbital elements, then converts to
% vectors in heliocentric frame.
[rxyzE0,vxyzE0]=OEtoXYZ(oeE(1),oeE(2),oeE(3),oeE(4),oeE(5),oeE(6),1);
[rxyzM0,vxyzM0]=OEtoXYZ(oeM(1),oeM(2),oeM(3),oeM(4),oeM(5),oeM(6),1);
[rxyzJ0,vxyzJ0]=OEtoXYZ(oeJ(1),oeJ(2),oeJ(3),oeJ(4),oeJ(5),oeJ(6),1);

%% Task 1
% Obtains rxyz and vxyz vectors for the planets on Mar 14, 2021 at 0159 UTC
% as well as the true anomalies
t1=[2021,3,14,01,59,00]; %specified time as a vector
ToF1=ToFfun(t1); %finding time of flight in AU
[rxyzE1,vxyzE1]=uToF(rxyzE0,vxyzE0,ToF1,1); %obtaining final vectors for Earth
[thE]=Tanomaly(rxyzE1,vxyzE1,1); %Earth true anomaly
[rxyzM1,vxyzM1]=uToF(rxyzM0,vxyzM0,ToF1,1); %obtaining final vectors for Mars
[thM]=Tanomaly(rxyzM1,vxyzM1,1);
[rxyzJ1,vxyzJ1]=uToF(rxyzJ0,vxyzJ0,ToF1,1); %obtaining final vectors for Jupiter
[thJ]=Tanomaly(rxyzJ1,vxyzJ1,1);
fprintf("Earth anomaly: %5.2f Mars anomaly %5.2f",thE,thM); %outputting true anomalies
fprintf("\n should be 70.87 and 291.3\n"); %reference values provided by Xu

%% Functions
% Organized here for ease of editing
function [rxyz,vxyz] = OEtoXYZ(a,e,i,O,w,th,mu)
%takes orbital elements and determines heliocentric-ecliptic vectors
p=a*(1-e^2);
rm=(p)/(1+e*cosd(th));
rpkw=[rm*cosd(th);rm*sind(th);0]; %perifocal vectors
vpkw=sqrt(mu/p)*[-sind(th);e+cosd(th);0];
R=[cosd(O)*cosd(w)-sind(O)*sind(w)*cosd(i),-cosd(O)*sind(w)-sind(O)*cosd(w)*cosd(i),sind(O)*sind(i);sind(O)*cosd(w)+cosd(O)*sind(w)*cosd(i),-sind(O)*sind(w)+cosd(O)*cosd(w)*cosd(i),-cosd(O)*cosd(i);sind(w)*sind(i),cosd(w)*sind(i),cosd(i)];
rxyz=R*rpkw; %heliocentric-ecliptic vectors
vxyz=R*vpkw;
end

function [r1,v1] = uToF(r0,v0,ToF,mu)
% Universal ToF using z (for all orbit types)
r0m=norm(r0); %radius magnitude
v0m=norm(v0); %velocity magnitude
ene=v0m^2/2-mu/r0m; %orbital energy
a=-mu/(2*ene); % semi-major axis
xo=ToF*sqrt(mu)/a; %initial guess for x
zo=xo^2/a; %initial guess for z
% Selecting anonymous functions used during iteration based on orbit
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
% iterative solver
while 1
    t=((xo^3*S(zo))+(dot(r0,v0)/sqrt(mu)*xo^2*C(zo))+(r0m*xo*(1-zo*S(zo))))/sqrt(mu);
    rn=xo^2*C(zo)+dot(r0,v0)/sqrt(mu)*xo*(1-zo*S(zo))+r0m*(1-zo*C(zo));
    if abs((ToF-t)/ToF) <= 1e-5 %terminating if difference between times is small enough
        x=xo;
        break
    end
    xo=xo+(ToF-t)/(rn/sqrt(mu)); % newton-step to next x guess
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

% function [ToF] = TOFcalc(t1) %might just use anonymous function
% %takes datetime vector input and finds ToF in heliocentric TU
% t0=datetime(2000,1,1,11,58,0); %setting initial time to the J2000 parameter
% ToF=etime(timevec(t1),timevec(t0))/5.0226757e6; %should remove datetime and timevec for actual use
% end

function [th] = Tanomaly(rxyz,vxyz,mu)
% takes heliocentric-ecliptic vectors and outputs the true anomaly in deg
rm=norm(rxyz);
ev=1/mu*((dot(vxyz,vxyz)-mu/rm)*rxyz-(dot(rxyz,vxyz)*vxyz));
e=norm(ev);
th=acosd(dot(ev,rxyz)/(e*rm));
if dot(rxyz,vxyz)<0
    th=360-th;
end
end

%% Gauss Orbit
% inputs r0, ToF, mu, r1, zg
function [v0,v1] = Gorb(r0,r1,ToF,mu,zg)
r0m=norm(r0); %initial radius magnitude
r1m=norm(r1); %final radius magnitude
Dths=acosd(dot(r0,r1)/(r0m*r1m)); %delta theta (short path), in degrees

As=sqrt(r0m*r1m)*sind(Dths)/sqrt(1-cosd(Dths)); %short path A
% Selecting anonymous functions used during iteration based on orbit
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
dfs=(-sqrt(mu)*xs/(r0m*r1m))*(1-zs*S(zs));
dgs=1-xs^2/r1m*C(zs);
v0s=(r1-fs*r0)/gs; %outputting short path initial and final velocity vectors
v1s=(dgs*r1-r0)/gs;
end