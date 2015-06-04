function [deltaV_hohman , Hohman_time] = Hohman_Transfer(R1, R2, mu)
%Takes two radius's assuming transfer between two circular orbits
%that are on the same plane (so delta_i=delta_omega=0) Also, takes the
%mu parameter depending on the body

a = 0.5*(R1 + R2); %semi-major axis for transfer ellipse
v_c1 = sqrt(mu/R1); %Circular Orbital Velocity1
v_c2 = sqrt(mu/R2); %Circular Orbital Velocity2
roe = R2/R1; %defined ratio
if roe > 11.9
    disp('Note: R2/R1>11.9; Bi-Elliptic Transfer could be more optimal')
end
%From Lecture Notes:
deltaV_1 = v_c1*(sqrt(2*roe/(1+roe))-1);
deltaV_2 = v_c1*(sqrt(1/roe)-sqrt(2/(roe*(1+roe))));


deltaV_hohman = deltaV_1 + deltaV_2;
Hohman_time = (a)^(3/2)*pi/sqrt(mu);

end 