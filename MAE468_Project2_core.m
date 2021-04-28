% This is the core code file for the MAE468 Project 2 submission
% The team consists of Joseph Barragree, Sarah Polickoski, Micajah
% Schweikert, and Stephen Ward.

%% Notes
%source used for some constants: http://www.dept.aoe.vt.edu/~lutze/AOE2104/consts.pdf
%source for Mars and Earth constants https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
%ADCS reaction wheel source, using HR 12-25 https://satcatalog.com/datasheet/Honeywell%20-%20HR12-25.pdf
%ADCS thruster source, using MR-103J https://www.rocket.com/sites/default/files/documents/In-Space%20Data%20Sheets%204.8.20.pdf
%Launch Vehicle source, using Falcon Heavy https://www.spacelaunchreport.com/falconH.html
%Launch Vehicle cost, Falcon Heavy https://www.faa.gov/about/office_org/headquarters_offices/ast/media/2018_ast_compendium.pdf

%% Housekeeping
% Run to remove figures, workspace variables and command window content
format compact
close all
clear
clc

%% Orbital Variable Initialization
% Includes orbital elements, locations, constants, anonymous functions, etc
oeE=[1.000000,0.01671,0.00005,-11.26064,114.20783,-2.48284]; %Earth
oeM=[1.523662,0.093412,1.85061,49.57854,286.4623,19.41248]; %Mars
% a,e,i,O,w,th is naming convention used in functions
t0=datetime(2000,1,1,11,58,0); %setting initial time to the J2000 parameter
ToFfun=@(tt) etime(datevec(tt),datevec(t0))/5.0226757e6; %anonymous function for finding ToF in AU
zg=18; %initial guess for z (for Gauss Orbit). Should be within the range of +-(2pi)^2 
% NOTE: Use 18 if elliptical, 0 if parabolic, -18 if hyperbolic
vCAtoM=@(v) v*149597870.7/5.0226757e6; %km/s (AU/TU sun to km/s)
muM=4.2828e4; %Mars mu, km^3/s^2
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
fprintf("Earth at Departure\n\t Position: %5.4f, %5.4f, %5.4f AU\n\t Velocity: %5.4f, %5.4f, %5.4f AU/TU\n",rxyzE1(1),rxyzE1(2),rxyzE1(3),vxyzE1(1),vxyzE1(2),vxyzE1(3));
fprintf("Mars at Arrival\n\t Position: %5.4f, %5.4f, %5.4f AU\n\t Velocity: %5.4f, %5.4f, %5.4f AU/TU\n",rxyzM2(1),rxyzM2(2),rxyzM2(3),vxyzM2(1),vxyzM2(2),vxyzM2(3));

%% Task c
% Spacecraft interplanetary transfer parameters
fprintf("\n---TASK c---\n");
[vxyz1,vxyz2]=Gorb(rxyzE1,rxyzM2,ToFfun(t0+days(190)),1,zg); %spacecraft heliocentric departure and arrival velocities
[oeSC1(1),oeSC1(2),oeSC1(3),oeSC1(4),oeSC1(5),oeSC1(6)]=XYZtoOE(rxyzE1,vxyz1,1); %Spacecraft departure orbital elements
[~,~,~,~,~,oeSC2]=XYZtoOE(rxyzM2,vxyz2,1); %Spacecraft arrival true anomaly
fprintf("Spacecraft RoI Orbital Elements\n\t Eccentricity: %5.3f\n\t Semi-major axis: %5.3f AU\n\t Inclination: %5.3f degrees\n\t RAAN: %5.3f degrees\n\t Argument of periapsis: %5.3f degrees\n\t Departure true anomaly: %5.3f degrees\n\t Arrival true anomaly: %5.3f degrees\n",oeSC1(1),oeSC1(2),oeSC1(3),oeSC1(4),oeSC1(5),oeSC1(6),oeSC2);
fprintf("Spacecraft at Earth RoI departure\n\t Heliocentric Velocity: %5.4f, %5.4f, %5.4f AU/TU\n",vxyz1(1),vxyz1(2),vxyz1(3));
fprintf("Spacecraft at Mars RoI Arrival\n\t Heliocentric Velocity: %5.4f, %5.4f, %5.4f AU/TU\n",vxyz2(1),vxyz2(2),vxyz2(3));

%% Task d
% Find delta v and escape velocity for the trajectory assuming Earth
% launch and 400km circular Mars parking orbit
fprintf("\n---TASK d---\n");
alt=400; %spacecraft parking orbit altitude in km
vinf1=vCAtoM(norm(vxyz1-vxyzE1)); %finding vinf at Earth in km/s
vinf2=vCAtoM(norm(vxyz2-vxyzM2)); %finding vinf at Mars in km/s
dvM=sqrt(muM/(3389.5+alt))-sqrt(2*(vinf2^2/2+muM/(3389.5+alt)));%finding dv for Mars 400km altitude parking orbit
fprintf("Mars vinf to parking orbit dv: %5.4f km/s\n",abs(dvM)); %displaying delta v to get Mars parking orbit

%% Transmitter sizing
% Size transmitter based on longest distance, data transfer with 2kbps
% additional engineering downlink, 10^-6 BER conv coding (Eb/N0=4.7), 10dB
% margin, efficiency based off of HW9 MGS
fprintf("\n---TASK e---\n");
CasIns=[11.83,14.46,32;55.9,57.8,365;3.1,3,3.6;27.7,9.25,1.5]; %Cassini Instrument parameters [W,kg,max kb/s]
DSN=[70,0.7,21]; %DSN [diameter,efficiency,noise temperature]
comm=[13.8e9,100e6]; %communications [frequency, bandwidth] in Hz
Mdist=[401e9,3389.5,227.923e6]; %Mars distance information [max dist from Earth m, planet radius m, average orbit distance km]
xmitt=[2.5,0.55,0,0]; %initializing transmitter information. [diameter m, efficiency,mass kg,power W]
Gain=@(D,f,n) -159.6+10*log10(n)+20*log10(D)+20*log10(f); %antenna gain in dB
Tloss=@(x,f) -147.6+20*log10(x)+20*log10(f)+1;%total path and other loss in dB

xmitt(3)=2.89*xmitt(1)^2+6.11*xmitt(1)-2.59; %calculating mass of transmitter in kg
xmitt(4)=10^((14.7-Gain(xmitt(1),comm(1),xmitt(2))-Gain(DSN(1),comm(1),DSN(2))+Tloss(Mdist(1),comm(1))+10*log10(1e3*sum(CasIns(:,3))+2000)-228.6+10*log10(DSN(3)))/10); %transmitter power in W
fprintf("Transmitter parameters\n\t Diameter: %3.1f m\n\t Efficiency: %3.2f\n\t Max power consumption: %4.2f W\n\t Mass: %4.2f kg\n",xmitt(1),xmitt(2),xmitt(4),xmitt(3));

%% Power system sizing
% Use calculated transmitter power consumption with Cassini instrument
% consumption and plug into planetary orbiter total power equation with 30% subsystem margin to determine
% solar panel sizing (assume Mars circular orbit for now). Also account for system degradation of 17yr mission.
% Find battery size accounting for degradation and appropriate power cycle life
[TlM,TnM]=Ltime(Mdist(2),alt,muM,0); %assuming th0=0 for most conservative estimate [time in light,time in night] in seconds
Pln=[sum(CasIns(:,1))+xmitt(4),sum(CasIns(:,1))]; %worst case payload power summation [power in light,power in night] both in W
Pln=Pln+1.3*((332.93*log(max(Pln))-1046.6)-max(Pln)); %total orbiter power consumption [power in light,power in night] both in W

Xln=[.8,.6]; %peak power transfer efficiencies light and night
PsaM=((Pln(1)*TlM)/Xln(1)+(Pln(2)*TnM)/Xln(2))/TlM; %solar array required power
PeolM=301*(149.596e6/Mdist(3))^2*0.77*cosd(0)*(1-0.005)^17; %end of life unit area power [W/m^2] for multijunction over 17 years
SolarArray=[PsaM/PeolM,0,0]; %initializing solar array information [array area m,array mass kg, battery mass kg]
SolarArray(2)=SolarArray(1)*4.0; %solar array mass in kg based on required area and rigid fold-out panels
SolarArray(3)=((Pln(2)*TnM)/(0.40*.97))/(55*60*60); %determining required battery mass based on night power consumption, 40%DoD, 97% transfer efficiency, and SED of 55W-hr/kg
fprintf("Power System Parameters\n\t Panel area: %4.2f m\n\t Array mass: %4.2f kg\n\t Battery mass: %4.2f kg\n",SolarArray(1),SolarArray(2),SolarArray(3));

%% Thermal system sizing
% Size thermal system based off of solar panels, arbitrary reasonable spacecraft size,
% and Mars circular orbit. Figure out radiator size if needed and heater
% power consumption. Note, transmitter will only run during the day
sfb=5.670374e-8; %stefan-boltzmann constant W/(m^2*K^4)
qMir=[162,120]*Mdist(2)^2/(Mdist(2)+alt)^2; %Mars max and min IR for parking orbit
GsM=1367*(149.596e6/Mdist(3))^2; %Solar at Mars based off of reference distances
Qwaste=[(sum(CasIns(:,1))+xmitt(4))/2+80,80]; %electronics waste heat [max,min] both in W
albM=0.29; %Mars albedo
Tlim=[30,10]+273; %spacecraft temperature limits [max,min] both in K <see paper for assumptions>
SCdim=[1.55,0,0]; %initializing spacecraft dimensions [radius m,area in m^2,area out m^2]
SCdim=[SCdim(1),pi()*SCdim(1)^2,4*pi()*SCdim(1)^2]; %calculating the areas
ae=[0.30,0.03,0.8,0.78,0.17,0.92]; %surface radiation factors [SCgold abs,SCgold emis,solar panel abs,radiator emis,SCpaint abs,SCpaint emis] 
%using vapor deposited gold SC coating and 5mil alumized teflon radiators with some of the spacecraft painted white with Z93 paint
cfrac=0.216; %fraction of spacecraft surface that is painted white instead of gold coated

Qsolalb=(1+albM)*GsM*((ae(5)*cfrac+ae(1)*(1-cfrac))*SCdim(2)+ae(3)*SolarArray(1)); %total solar Qin including albedo and solar panels [max,min] in W
Qir=qMir*(SCdim(2)+SolarArray(1)); %planet IR Qin [max,min] in W
%Temps=[((Qsolalb+max(Qir)+max(Qwaste))/(SCdim(3)*ae(2)*sfb))^(1/4),((min(Qir)+min(Qwaste))/(SCdim(3)*ae(2)*sfb))^(1/4)]; %SC temperatures without radiator [max,min] both in K
Arad=[(Qsolalb+max(Qir)+max(Qwaste)-SCdim(3)*(ae(6)*cfrac+ae(2)*(1-cfrac))*sfb*max(Tlim)^4)/(ae(4)*sfb*max(Tlim)^4),(min(Qir)+min(Qwaste)-SCdim(3)*(ae(6)*cfrac+ae(2)*(1-cfrac))*sfb*min(Tlim)^4)/(ae(4)*sfb*min(Tlim)^4)]; %deployed radiator sizing for hot and cold cases
mTsys=[(0.24*cfrac+19.3e3*100e-9*(1-cfrac))*SCdim(3),max(Arad)*2.707]; %thermal system masses [mass of paint and gold coating kg,mass of radiators kg]
fprintf("Thermal System Parameters\n\t Deployed radiator area: %4.2f m^2\n\t Retracted radiator area: %4.2f m^2\n\t Thermal system mass: %4.2f kg\n",Arad(1),Arad(2),sum(mTsys));

%% ADCS Sizing
% Determine overall mass for craft based on all other components, add 60kg for CDS. Figure
% out approximate mass moments of inertia based on spherical (and homogeneous) craft, and
% deployed solar panels (using their mass). Determine grav-gradient and solar disturbances and size
% reaction wheels for a reasonable number of orbits maintaining within 5deg, then size thrusters.
% Figure out propellant mass for mission life with thrusters, then add margin.
SCmI=[sum(CasIns(:,2))+xmitt(3)+SolarArray(3)+mTsys(1)+60,0,0,0,0]; %initializing spacecraft sphere mass and mass moments of inertia [sphere mass,xmoment,ymoment,zmoment,total mass]
SCmI(2)=2/5*SCmI(1)*SCdim(1)^2+((SCdim(1)+0.25+SolarArray(1)/8)^2+((SolarArray(1)/4)^2+0.05^2)/12)*SolarArray(2)+((SCdim(1)-0.25)^2+(3^2+0.5^2)/12)*mTsys(2); %x axis moment
SCmI(3)=2/5*SCmI(1)*SCdim(1)^2+(0.05^2+2^2)/12*SolarArray(2)+((SCdim(1)-0.25)^2+(0.5^2+(max(Arad)*0.05)^2)/12)*mTsys(2); %y axis moment, see paper for assumptions
SCmI(4)=2/5*SCmI(1)*SCdim(1)^2+((SCdim(1)+0.25+SolarArray(1)/8)^2+((SolarArray(1)/4)^2+2^2)/12)*SolarArray(2)+(3^2+(max(Arad)*0.05)^2)/12*mTsys(2); %z axis moment

orbits=100; %variable to allow simple changing of orbits per momentum dump
SCdstrb=Tgrav(SCmI(4),SCmI(3),Mdist(2)+alt,muM,5);%symmetric solar panels, so only grav-gradient torque for 5deg attitude
SCdstrb(2)=SCdstrb(1)*(TlM+TnM)/4*0.707; %orbit total disturbance in N*m*s
[ADCS(1),ADCS(2)]=Mdump(orbits*SCdstrb(2),SCdim(1),10,202); %dumping stored momentum, some manual iteration required to find parameters to keep required thrust under 1N
ADCS(3)=4*(6)+8*(0.37)+(ADCS(2)/orbits*((60*60*24*365.2425*17)/(TlM+TnM))); %ADCS mass, accounts for 4 reaction wheels, 8 thrusters, and all necessary fuel
fprintf("ADCS System Parameters\n\t Stored momentum for %2.0f orbits: %5.2f N*m*s\n\t Thrust for 10 second momentum dump: %4.2f N\n\t Fuel required for mission life: %4.2f kg\n\t ADCS system mass: %5.2f kg\n",orbits,orbits*SCdstrb(2),ADCS(1),ADCS(2)/orbits*((60*60*24*365.2425*17)/(TlM+TnM)),ADCS(3));

%% Propulsion sizing
% Select engine parameters for an efficient but powerful enough chemical
% engine for Mars arrival burn. Based on the orbiter size, determine
% inert mass and fuel required.
propP=[exp(abs(dvM*1000)/(9.80665*223)),0.3]; %propulsion parameters. [MR,IMF] finding orbiter worst case MR based on Mars injection orbit with MR-104H thrusters
SCmI(5)=SCmI(1)+mTsys(2)+SolarArray(2)+ADCS(3); %total spacecraft mass
SCmI(6)=SCmI(5)*(propP(1)-1)*(1-propP(2))/(1-propP(1)*propP(2)); %orbiter propellant mass
SCmI(7)=propP(2)/(1-propP(2))*SCmI(6); %orbiter inert mass
SCmI(8)=sum(SCmI(6:7)); %orbiter total mass
fprintf("Propulsion System Parameters\n\t Payload mass: %6.2f kg\n\t Propellant mass: %6.2f kg \n\t Inert mass: %5.2f kg\n\t Total spacecraft mass: %6.2f kg\n",SCmI(5),SCmI(6),SCmI(7),SCmI(8));

%% Launch craft sizing
% Once mass is determined, investigate commercial launch options. If our
% craft is too heavy, figure out the fuel and mass for a two or three stage
% to orbit rocket.
fprintf("Launch Craft Selection\n\t Based on spacecraft mass and desired transfer orbit, a Falcon Heavy will\n\t be used to obtain the required Earth hyperbolic Vinf of %5.3f km/s\n",vinf1);

%% Price of mission
% Determine price of mission based on generic rates or specific
% component information.
cost=0.2*SCmI(8)+20*(17+190/365.2425)+270; %cost in M of dollars
fprintf("Total mission cost: $%4.3f Billion\n",cost/1000);

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

function [Tl,Tn] = Ltime(R,H,mu,th0)
% time in light and night for solar panel calculations
% inputs: R=planet radius [km] H=orbit altitude [km], 
% mu=planet mu [km^3/s^2], th0=eclipse angle (max use 0)
Lfrac=(pi()+2*asin((R*tan(acos(R./(R+H))))./((R+H)*cosd(th0))))./(2*pi()); %fraction of time in light
Tl=Lfrac.*(2*(R+H).^(3/2)*pi())/sqrt(mu); %time in light
Tn=(1-Lfrac).*(2*(R+H).^(3/2)*pi())/sqrt(mu); %time in night
end

function [Tg] = Tgrav(Iz,Iy,R,mu,th)
Tg=3/2*mu/R^3*abs(Iz-Iy)*sind(2*th);
end

function [Ft,mp] = Mdump(h,L,t,isp)
% Thruster sizing and propellant use for momentum dump
Ft=h/(L*t); %thruster force
mp=Ft*t/(isp*9.80665); %propellant mass
end
