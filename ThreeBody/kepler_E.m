function E = kepler_E(e, M) 
error = 1.e-8; %Set an error tolerance
%Select a starting value for E:
if M < pi 
    E = M + e/2; 
else
    E = M - e/2;
end
%Iterate on Equation until E is determined to within the error tolerance:
ratio = 1; 
while abs(ratio) > error 
    ratio = (E - e*sin(E) - M)/(1 - e*cos(E)); 
    E = E - ratio; 
end
end 
