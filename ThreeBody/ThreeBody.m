%-----------------------------------------------------------------------
%{
 This program uses kepler's equation for orbital bodies and newton's method
 to illustrate a spacecraft in an inclined motion around the moon. Also
 calculates the max velocity and acceleration experienced as a restricted
 three-body problem relevant to an interial rotating frame. 

 days -             converts days to seconds
 G -                universal graviational constant (km^3/kg/s^2)
 rmoon -            radius of the moon (km)
 rearth -           radius of the earth (km)
 r12 -              distance from center of earth to center of moon (km)
 m1,m2 -            masses of the earth and of the moon, respectively (kg)
 M -                total mass of the restricted 3-body system (kg)
 mu -               gravitational parameter of earth-moon system (km^3/s^2)
 mu1,mu2 -          gravitational parameters of the earth and of the moon,
                    respectively (km^3/s^2)
 Period_e,_m -      earth's and moon's periods of rotation (s)
 t_moon -           time for moon to orbit the earth (s)
 w_e,w_m -          angular velocity of earth and moon (rad/s)
 pi_1,pi_2 -        ratios of the earth mass and the moon mass, respectively,
                    to the total earth-moon mass
 W -                angular velocity of moon around the earth (rad/s)
 x1,x2 -            x-coordinates of the earth and of the moon, respectively,
                    relative to the earth-moon barycenter (km)

 User M-functions required: kepler_E
 User subfunctions required: none
%}
%-----------------------------------------------------------------------
warning('off','all');
clear all; close all; clc 
%Inputted Values--------------------------------------------------------
e = .65;%Eccentricity
a = 6542;%Semi-major Axis Length
i = [pi/2;pi/2;pi/2;pi/2] + pi/180*[60;60;60;60];%Inclination
%Standard values (obtained from JPL) and Calculations-------------------
G = 6.673e-20; 
rmoon = 1737.4; 
rearth = 6378.1; 
r12 = 384403; 
m1 = 5974e21; 
m2 = 7348e19; 
Period_e = 86400;               
Period_m = 2.3605776e6;    
t_moon = 2.3605776e6;  
w_e = 2*pi/Period_e;   
w_m = 2*pi/Period_m;    
M = m1 + m2;
pi_1 = m1/M; 
pi_2 = m2/M; 
mu1 = 398600.4418; 
mu2 = 4903.027779; 
mu = mu1 + mu2; 
W = sqrt(mu/r12^3); 
x1 = -pi_2*r12; 
x2 = pi_1*r12; 
b = a*sqrt(1-e^2);
%Tao = time delay from initial starting point---------------------------
tao=[0;23740.750;100;1000;200];
%Accuracy of Data(do not change unless experienced)---------------------
tend = Period_m;
span = 30;
t_span = 0:span:tend;  
Period = 2*pi*sqrt(a^3/mu2);
SALength = Period*3/span;
%Intializizing constraints----------------------------------------------
minspeed = 1000000000000000;
maxspeed = 0;
minaccel = 1000000000000000;
maxaccel = 0;

for n=1:size(t_span(:))
    t=t_span(n);
    time(n) = t;
    %Rotation Matrix of Moon--------------------------------------------
    R=[cos(w_m*t) sin(w_m*t) 0; -sin(w_m*t) cos(w_m*t) 0; 0 0 1];
    Rdot=[-w_m*sin(w_m*t) w_m*cos(w_m*t) 0;-w_m*cos(w_m*t) -w_m*sin(w_m*t) 0; 0 0 0];
    %Position, Velocity of Moon in ECI----------------------------------
    Moon_ECI=[r12*cos(w_m*t);r12*sin(w_m*t);0];
    Moon_V_ECI=[-r12*w_m*sin(w_m*t); r12*w_m*cos(w_m*t);0];
    Speed_ECI(n)=norm(Moon_V_ECI);
    x_ECI(n)=Moon_ECI(1);
    y_ECI(n)=Moon_ECI(2);
    z_ECI(n)=Moon_ECI(3);
    %-------------------------------------------------------------------
    M = sqrt(mu2/a^3)*(t-tao(1));
    E = kepler_E(e,M); 
    thetha=2*atan(sqrt((1+e)/(1-e))*tan(E/2));
    displayt = thetha*180/pi;
    theR = a*(1-e^2)/(1+e*cos(thetha));
    x_1 = theR*cos(thetha);
    y_1 = theR*sin(thetha);
    z_1 = 0;
    r_1=[x_1;y_1;z_1];
    tt = 90*pi/180;
    %Conversion to Moon Centered/Three Body-----------------------------
    Rot2 = [cos(i(1)) 0 -sin(i(1)); 0 1 0 ; sin(i(1)) 0 cos(i(1))];
    Rot3 = [cos(tt) sin(tt) 0; -sin(tt) cos(tt) 0; 0 0 1];
    r_F = Rot2*[x_1;y_1;z_1];
    x = pi_1*r12 - r_F(1);
    y = r_F(2);
    z = r_F(3);
    rs = [x;y;z];
    r1 = [x+pi_2*r12;y;z];
    r2 = [x-pi_1*r12;y;z];
    r_ECI = Moon_ECI +R*r_F;
    xSC_ECI(n)=r_ECI(1);
    ySC_ECI(n)=r_ECI(2);
    zSC_ECI(n)=r_ECI(3);
    a_sc = (-mu1/(norm(r1))^3)*r1 - (mu2/(norm(r2))^3)*r2;
    Accel_SC(n)=norm(a_sc);
    %Min/Max Acceleration-----------------------------------------------
    if Accel_SC(n) < minaccel
        minaccel = Accel_SC(n);
    end
    if Accel_SC(n) > maxaccel
        maxaccel = Accel_SC(n);
    end
    %Velocity in MCI---------------------------------------------------
    speed = sqrt(mu2*(2/norm(r_1)-1/a));
    Speed_SC(n) = speed; 
    fpa = atan((e*sin(thetha))/(1+e*cos(thetha)));
    runit = r_F/norm(r_F);
    runitp1 = Rot3*(r_1/norm(r_1));
    runitp = Rot2*runitp1;
    vper = speed*cos(fpa);
    vr = speed*sin(fpa);
    v_SC = vr*runit+vper*runitp;
    %Common Oribital Values---------------------------------------------
    h = cross(r_ECI,v_SC);
    khat = [0;0;1];
    ni = [-h(2);h(1);0];
    if ni(2)>0
        RAAN = 180/pi*acos(ni(1)/norm(ni));
    end
    if ni(2)<0 
        RAAN = 180/pi*(2*pi - acos(ni(1)/norm(ni)));
    end
    ECC = cross(v_SC,h)/mu2 - r_ECI/norm(r_ECI);
    AOP = acos(dot(ni,ECC)/(norm(ni)*norm(ECC)))*180/pi;
    %Max/Min Speed------------------------------------------------------
    if Speed_SC(n) < minspeed
        minspeed = Speed_SC(n);
    end
    if Speed_SC(n) > maxspeed
        maxspeed = Speed_SC(n);
    end
    
    %Euler Angles-------------------------------------------------------
    rsun = [1;0;0];
    r1a = [x_1;y_1;0];
    rss = [r_F(1); 0 ; r_F(3)];
    rs_1 = [r_F(1); 0 ; r_F(3)] + [0;0;-rmoon];
    if y_1 == 0 || y_1 > 0
        Aangle1 = acos(dot(r1a,rsun)/norm(r1a)) + pi;
    end
    if y_1 < 0
        Aangle1 = pi - acos(dot(r1a,rsun)/norm(r1a));
    end
    Aangle3 = acos(dot(rsun,rs_1)/norm(rs_1))/2;
    Aangle4 = acos(dot(rss,rs_1)/(norm(rs_1)*norm(rss)));
    Aangle2 = pi/2 + Aangle3 + Aangle4;
    if Aangle2 > 2.4
         Aangle2 = 2.4;
    end
    R3_Aangle1 = [cos(Aangle1) sin(Aangle1) 0; -sin(Aangle1) cos(Aangle1) 0; 0 0 1]; 
    R2_Aangle2 = [cos(Aangle2) 0 -sin(Aangle2); 0 1 0 ; sin(Aangle2) 0 cos(Aangle2)];
    Q = R*(R2_Aangle2*(Rot2*R3_Aangle1));
    [phi thetha2 psi] = dcm_to_euler(Q);
    AA1(n) = Aangle1;
    AA2(n) = Aangle2;
    phii(n) = phi;
    thetha22(n) = thetha2;
    psii(n) = psi;
end
%-----------------------------------------------------------------------
figure(1);
r10 = rearth * ones(50, 50);
[th, phi] = meshgrid(linspace(0, 2*pi, 50), linspace(-pi, pi, 50));
[x,y,z] = sph2cart(th, phi, r10);
surf(peaks);
surface(x,y,z,'FaceColor', 'none')
hold on;
plot3(x_ECI,y_ECI,z_ECI,xSC_ECI,ySC_ECI,zSC_ECI);
set(gca, 'ylim', [-400000 400000], 'xlim', [-400000 400000], 'zlim', [-400000 400000]);
title('Moon and Spacecraft Trajectory in ECI frame(km)')
xlabel('x')
ylabel('y')
zlabel('z')
grid on
hold on;
%----------------------------------------------------------------------
figure(2);
plot(time(1:SALength), Speed_SC(1:SALength));
title('Speed of S/C vs. Time')
xlabel('Time(s)')
ylabel('Speed(km/s)')
%----------------------------------------------------------------------
figure(3);
plot(time(1:SALength), Accel_SC(1:SALength));
title('Acceleration of S/C vs. Time')
xlabel('Time(s)')
ylabel('Acceleration(km^2/s)')
%----------------------------------------------------------------------
figure(4);
plot(time(1:SALength), phii(1:SALength),time(1:SALength), thetha22(1:SALength),time(1:SALength), psii(1:SALength));
title('Euler Angles of S/C vs. Time')
xlabel('Time(s)')
ylabel('Phi,thetha,psi')
legend('phi','thetha','psi')
set(gca, 'ylim', [0 360])
%----------------------------------------------------------------------
figure(5);
plot(time(1:SALength), AA1(1:SALength), time(1:SALength), AA2(1:SALength))
fprintf('For the Spacecraft\n Minimum Speed is = %d\n Maximum Speed is = %d\n Minimum Acceleration is = %d\n Maximum Acceleration is = %d\n', minspeed, maxspeed, minaccel,maxspeed);