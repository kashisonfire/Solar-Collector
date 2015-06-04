% Calculating mass of fuel used from deltaV
function [mf,m0_new] = DeltaV_to_mfuel(deltaV, v_e, m0)
% Returns the mf fuel used and new mass from a cetain deltaV given v_e and 
% intial mass
% Note: Following derived from deltaV = -V_e*ln(m(t)/m0) and m(t)=m0-mf(t)
mf = m0*(1-exp(-deltaV/v_e));
m0_new = m0 - mf;
end