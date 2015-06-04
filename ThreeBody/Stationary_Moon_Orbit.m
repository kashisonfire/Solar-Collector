clear all; close all; clc 

days = 24*3600; 
G = 6.6742e-20; 
rmoon = 1737; 
rearth = 6378; 
r12 = 384400; 
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

mu1 = 398600; 
mu2 = 4903.02; 
mu = mu1 + mu2; 

W = sqrt(mu/r12^3); 
x1 = -pi_2*r12; 
x2 = pi_1*r12; 
tao=[0;0;0;0];

e = [.576;.65];
a = [5400.3;6542];
i = [pi/2;pi/2;pi/2;pi/2] + pi/180*[60;60;180;270];
t_span = 0:30:Period_m;  

for n=1:size(t_span(:))
    t=t_span(n);
    M1 = sqrt(mu2/a(1)^3)*t;
    E1 = kepler_E(e(1),M1); 
    thetha1=2*atan(sqrt((1+e(1))/(1-e(1)))*tan(E1/2));
    theR1 = a(1)*(1-e(1)^2)/(1+e(1)*cos(thetha1));
    x1 = theR1*cos(thetha1);
    y1 = theR1*sin(thetha1);
    z1 = 0;

    Rot1 = [cos(i(1)) 0 -sin(i(1)); 0 1 0 ; sin(i(1)) 0 cos(i(1))];
    r1 = Rot1*[x1;y1;z1]; 
    xSC1(n)=r1(1);
    ySC1(n)=r1(2);
    zSC1(n)=r1(3);   
    
    M2 = sqrt(mu2/a(2)^3)*t;
    E2 = kepler_E(e(2),M2); 
    thetha2=2*atan(sqrt((1+e(2))/(1-e(2)))*tan(E2/2));
    theR2 = a(2)*(1-e(2)^2)/(1+e(2)*cos(thetha2));
    x2 = theR2*cos(thetha2);
    y2 = theR2*sin(thetha2);
    z2 = 0;
    
    Rot2 = [cos(i(1)) 0 -sin(i(1)); 0 1 0 ; sin(i(1)) 0 cos(i(1))];
    r2 = Rot2*[x2;y2;z2]; 
    xSC2(n)=r2(1);
    ySC2(n)=r2(2);
    zSC2(n)=r2(3);
    
%     Rot3 = [cos(i(3)) 0 -sin(i(3)); 0 1 0 ; sin(i(3)) 0 cos(i(3))];
%     r3 = Rot3*[x;y;z]; 
%     xSC3(n)=r3(1);
%     ySC3(n)=r3(2);
%     zSC3(n)=r3(3);
%     
%     Rot4 = [cos(i(4)) 0 -sin(i(4)); 0 1 0 ; sin(i(4)) 0 cos(i(4))];
%     r4 = Rot4*[x;y;z]; 
%     xSC4(n)=r4(1);
%     ySC4(n)=r4(2);
%     zSC4(n)=r4(3);
    
end
r10 = rmoon * ones(50, 50); 
[th, phi] = meshgrid(linspace(0, 2*pi, 50), linspace(-pi, pi, 50));
[x,y,z] = sph2cart(th, phi, r10);
surf(peaks);
surface(x,y,z,'FaceColor', 'none')
hold on;
plot3(xSC1,ySC1,zSC1,xSC2,ySC2,zSC2);%xSC3,ySC3,zSC3,xSC4,ySC4,zSC4);
set(gca, 'ylim', [-10000 10000], 'xlim', [-10000 10000], 'zlim', [-10000 10000]);
hold off;
Grid on;