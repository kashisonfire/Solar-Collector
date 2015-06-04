%% Main file: integration + plotting

%% Preamble
% Some constants: units: rad, sec, m, kg
d2r = 180.0/pi;
r2d = pi/180.0;

minute = 60;

Mass = 1000;
rho = 2810;
Radius = .5;
height = Mass/(pi*Radius^2);


Inertia = [1/12*Mass*height^2+1/4*Mass*Radius^2 0.0 0.0; 
           0.0 1/12*Mass*height^2+1/4*Mass*Radius^2 0.0; 
           0.0 0.0 1/2*Mass*Radius^2];
       
%% Initial conditions and integration
i = 60*r2d;
bm = [cos(i) 0 -sin(i); 0 1 0; sin(i) 0 cos(i)];

b10 = [cos(i);0;-sin(i)];
b20 = [0;1;0];
b30 = [sin(i);0;cos(i)];
ang_vel0 = [0.01; 0.5; 0.01];
X0 = [b10; b20; b30; ang_vel0];

Tspan = [0 2.0*minute];

[T,Y] = ode45(@free_rigid_body,Tspan,X0);

b1 = Y(:,1:3);
b2 = Y(:,4:6);
b3 = Y(:,7:9);
ang_vel = Y(:,10:12);

%% Invariant ellipsoids

Npts = 100;
% Angular momentum ellipsoid
H = norm(Inertia*ang_vel0);
H1 = H/Inertia(1,1); % assumes that we are in principal axis.
H2 = H/Inertia(2,2);
H3 = H/Inertia(3,3);
[E1x, E1y, E1z] = ellipsoid(0.0, 0.0, 0.0, H1, H2, H3, Npts);

% Energy ellipsoid
K = 0.5*dot(ang_vel0, Inertia*ang_vel0);
K1 = sqrt(2.0*K/Inertia(1,1));
K2 = sqrt(2.0*K/Inertia(2,2));
K3 = sqrt(2.0*K/Inertia(3,3));
[E2x,E2y,E2z]= ellipsoid(0.0, 0.0, 0.0, K1, K2, K3, Npts);

%% Plotting

% Options
B = eye(3,3);
O = zeros(3,1);
S1 = r2d*max([H1, K1]); S2 = r2d*max([H2, K2]); S3 = r2d*max([H3, K3]); 
S = 1.1*max([S1, S2, S3]);
range = 1.3*[-S, S, -S, S, -S, S];
crange = [-S, S];

Ang_vel = ang_vel; % will be used to compute the ang_vel in inertial frame
e2x = E2x;
e2y = E2y;
e2z = E2z;

% Attitude movie
for i = 1:size(T),
    h=figure(2);
    clf;
    
    A = [b1(i,:);b2(i,:);b3(i,:)];
    Ai = A';
    hv = Inertia*(ang_vel(i,:)')*r2d; % factor r2d used for graphics scaling
    Hv = Ai*hv;
    
%---Body frame viewpoint     
    subplot(121)
    
    % Body frame
    h1 = quiver3( O, O, O, B(:,1), B(:,2), B(:,3), S);
    hold on;
    
    % Inertial frame 
    h2 = quiver3( O, O, O, Ai(:,1), Ai(:,2), Ai(:,3), S);
    
    % Polhode
    h3 = plot3( ang_vel(1:i,1)*r2d, ang_vel(1:i,2)*r2d, ang_vel(1:i,3)*r2d, '-r');
    set(h3,'LineWidth',2);
    
    % Angular velocity vector
    h4= quiver3( [0.0, 0.0],        [0.0, 0.0],        [0.0, 0.0],...
                 [0.0, ang_vel(i,1)]*r2d, [0.0, ang_vel(i,2)]*r2d,...
                 [0.0, ang_vel(i,3)]*r2d, 'g');
    set(h4,'LineWidth',2);
    
    % Angular momentum     
    h5 = quiver3( [0.0, 0.0],         [0.0, 0.0],         [0.0, 0.0],...
             [0.0, hv(1)], [0.0, hv(2)], [0.0, hv(3)], 'k');
    set(h5,'LineWidth',2)
         
    % Kinetic energy ellipsoid
    h6 =surf(r2d*E2x, r2d*E2y, r2d*E2z, -S/2.0+(S/2.0).*E2z);
    set(h6,'LineStyle','none','FaceAlpha',0.5);
    
    % Angular momentum ellipsoid
    h7 = surf(r2d*E1x, r2d*E1y, r2d*E1z, S/2.0+(S/2.0).*E1z);
    set(h7,'LineStyle','none','FaceAlpha',0.5);

    % Axis options
    grid on;
    axis equal;
    axis(range);
    caxis(crange);
    view(130, 30);
        
    % Labeling
    legend = ['body frame in blue, inertial frame in green\newline',...
              'angular velocity in green, angular momentum in black\newline',...
              'polhode in red, knietic energy ellipsoid in orange)'];
    title('Body frame viewpoint');
    text(-5*S,-5*S,legend);
    xlabel('b_1');
    ylabel('b_2');
    zlabel('b_3');
    hold off;
    
%---Inertial frame viewpoint 
    subplot(122)
        % Body frame
    h1 = quiver3( O, O, O, b1(i,:)', b2(i,:)', b3(i,:)', S);
    hold on;
    
    % Intertial frame 
    h2 = quiver3( O, O, O, B(:,1), B(:,2), B(:,3), S);
    
    % Polhode
    Ang_vel(i,:) = ang_vel(i,:)*A;
    h3 = plot3( Ang_vel(1:i,1)*r2d, Ang_vel(1:i,2)*r2d, Ang_vel(1:i,3)*r2d, '-r');
    set(h3,'LineWidth',2);
    
    % Angular velocity vector
    h4= quiver3( [0.0, 0.0],        [0.0, 0.0],        [0.0, 0.0],...
                 [0.0, Ang_vel(i,1)]*r2d, [0.0, Ang_vel(i,2)]*r2d,...
                 [0.0, Ang_vel(i,3)]*r2d, 'g');
    set(h4,'LineWidth',2);
    
    % Angular momentum     
    h5 = quiver3( [0.0, 0.0],         [0.0, 0.0],         [0.0, 0.0],...
             [0.0, Hv(1)], [0.0, Hv(2)], [0.0, Hv(3)], 'k');
    set(h5,'LineWidth',2)
         
    % Kinetic energy ellipsoid
    for j=1:Npts+1,
        for k=1:Npts+1
            E2 = [E2x(j,k); E2y(j,k); E2z(j,k)];
            e2 = Ai*E2;
            e2x(j,k) = e2(1);
            e2y(j,k) = e2(2);
            e2z(j,k) = e2(3);
        end;
    end;
    h6 =surf(r2d*e2x, r2d*e2y, r2d*e2z, -S/2.0+(S/2.0).*E2z);
    set(h6,'LineStyle','none','FaceAlpha',0.5);
    
    % Angular momentum ellipsoid
%    h7 = surf(r2d*E1x, r2d*E1y, r2d*E1z, S/2.0+(S/2.0).*E1z);
%    set(h7,'LineStyle','none','FaceAlpha',0.5);

    % Axis options
    grid on;
    axis equal;
    axis(range);
    caxis(crange);
    view(130, 30);
    
    % Labeling
    title('Inertial frame viewpoint');
    xlabel('x ');
    ylabel('y ');
    zlabel('z ');
    hold off;

end;
%movie2avi(F,'attitude_movie.avi')
   
