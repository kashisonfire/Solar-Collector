%% free rigid body attitude dynamics
% the state X is [b1, b2, b3, w], that is:
% 1) coords of the body frame unit vector in inertial space 
% 2) coords of the angular velocity in body frame 
% Units are rad, sec, m, kg
function dX = free_rigid_body(t,X)

% Initializations
dX=zeros(12,1);

% Inertia (!!! CAREFUL: match this to the defined Inertia in 'free_righid_body.m');
Mass = 1000;
rho = 2810;
Radius = .5;
height = Mass/(pi*Radius^2);

I1 = 1/12*Mass*height^2+1/4*Mass*Radius^2;
I2 = 1/12*Mass*height^2+1/4*Mass*Radius^2;
I3 = 1/2*Mass*Radius^2;

% Attiutde rotation
ang_vel = [X(10); X(11); X(12)]; % in body frame.
A = [X(1), X(2), X(3);...
     X(4), X(5), X(6);... 
     X(7), X(8), X(9)];          % direction cosine matrix

% angular velocity in inertial frame
OM =  A'*ang_vel;

% Attitude equations
dX(1:3) = cross(OM,X(1:3));
dX(4:6) = cross(OM,X(4:6));
dX(7:9) = cross(OM,X(7:9));

% Euler's equations
dX(10) = ((I2-I3)/I1)*X(11)*X(12);
dX(11) = ((I3-I1)/I2)*X(12)*X(10);
dX(12) = ((I1-I2)/I3)*X(10)*X(11);