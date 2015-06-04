%Testing Transfer Functions
clear all
clc
%Model parameters
G             = 6.673e-20;            %Universal gravitational constant in km^3/(kg*(s^2)
mu_earth      = 3.986004418e5;        %Gravitational parameter of Earth in km^3*s^-2
mu_moon       = 4.9027779e3;          %Gravitational parameter of the moon in km^3*s^-2
Period_E      = 86400;                %Earth's Period of Rotation in s
Period_M      = 2.3605776e6;          %Moon's Period of Rotation in s
t_moon        = 2.3605776e6;          %Time for moon to orbit Earth in s
w_E           = 2*pi/Period_E;        %Angular Velocity of Earth in rad/s
w_M           = 2*pi/Period_M;        %Angular Velocity of Moon in rad/s
r_E           = 6378.1;               %Radius of Earth in km
r_M           = 1737.4;               %Radius of Moon in km
d_m           = 384403;               %Center-Center distance btwn Earth and Moon in km
t_span        = 0:10:Period_M;        %Time Span in s


R_LEO         = 8378.1;               %Radius of Low Earth Orbit in km
R_GSO         = 42164;                %Apogee of Geosynchronous Transfer orbit and Radius of Geosynchronous orbits in km
i_M2E         = 5.145;                %Inclination (degrees) of Moon's Orbit to Earth
i_initial     = 27;                   %Where the Launch Vehicle Places us (27 for Falcon 9)
m0            = 4000;                 %Initial mass (4s/c) after separation from vehicle
I_sp          = 250;                  %Isp of rockets on s/c in seconds
g             = 9.81;                 %m/s^2
v_e           = I_sp*g;               %Definition of v_e from Isp (m/s)
v_e_km        = v_e/1000;             %v_e in km/s for one rocket
v_e_km4       = 4*v_e_km;             %Used for when four satellites are using their burns at once(4*Thrust=4*dm/dt*v_e) from F=ma for calculation of mass of fuel required

%%%%%%%
%%%%%%%
%Make eccentricity zero (circulat orbit) for hohman transfer (not sure if
%needed) but the launch vehicle puts us in GTO, which is ellipitical(e=.46)
v_cGSO = sqrt(mu_earth/R_GSO);         %velocity required for Geosynchronous orbit
v_aGTO = v_cGSO*sqrt(1-2*.46+.46^2);   %apogee velocity of GTO
deltaV_eto0 = v_cGSO-v_aGTO;
[mf_eto0 , m0] = DeltaV_to_mfuel(deltaV_eto0, v_e_km4, m0);
disp('DeltaV (km/s) for change to circular orbit is: ')
disp(deltaV_eto0)
disp('Mass of fuel used for change to circular orbit is: ')
disp(mf_eto0)
%%%%%%% Optional? Although Hohman Transfer is in between two circular orbits
%%%%%%%

%Changes inclination to the moons orbit inclination from our launch vehicle
%inclination
deltaV_inc = Incl_transfer(i_initial,i_M2E,R_GSO,mu_earth);
%Returns the deltaV needed for the inclination(or omega) change and also
%display a notice to the console if the change is greater than 38.94
%degrees, which is where other plane changes could be more optimal

%The Hohman Transfer from GSO to Moon
[deltaV_h, Transfer_t] = Hohman_Transfer(R_GSO,d_m,mu_earth); 
%Note: Returns two results delta V and transfer time. Also, displays a
%notice to the console if R2/R1 is bigger than 11.9 where other types of
%transfers could be more optimal


%Assuming that the 4 s/c's can travel together and fire all their rockets at
%once for the following manuevers (4*Thrust=4*-dm/dt*v_e)
[mf_inc , m01] = DeltaV_to_mfuel(deltaV_inc, v_e_km4, m0);
%Function returns the mass of fuel used and new mass from given deltaV and
%initial mass
[mf_hohman, m02] = DeltaV_to_mfuel(deltaV_h, v_e_km4, m01);

%Display results
disp('DeltaV (km/s) for inclination change near earth is: ')
disp(deltaV_inc)
disp('Total fuel (kg) used for inclination change near earth is: ')
disp(mf_inc)
disp('DeltaV (km/s) for hohman is: ')
disp(deltaV_h)
disp('Total fuel used (kg) for hohman is: ')
disp(mf_hohman)
disp('Transfer time (hours) for hohman is: ')
disp(Transfer_t/3600)


disp('Total fuel (kg) used with initial eccentricity change is: ')
disp(mf_hohman+mf_inc+mf_eto0)