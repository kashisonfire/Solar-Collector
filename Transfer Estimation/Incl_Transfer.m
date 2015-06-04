%INCLINATION/OMEGA IMPULSE CHANGE
function deltaV = Incl_Transfer(i1, i2, R, mu)
%i1, i2 in degrees, R in km, mu in km^3/s^2

%Takes a circular orbit and returns the deltaV for an inclination (or
%omega) change to another circular orbit with all other orbital elements
%remaining the same except the change in inclination (or omega)

v_1 = sqrt(mu/R); %Circular Orbital Velocity (km/s)
delta_i = abs(i2-i1); %Change in inclination (deg)
if delta_i > 38.94
    disp('NOTE: Angle Change difference is bigger than 38.94 degrees.')
    disp('May want to consider different plane change like bi elliptic for optimization')
end
delta_i_rad = delta_i*pi/180; %Change in inclination in rad

deltaV = 2*v_1*sin(delta_i_rad/2);
end 