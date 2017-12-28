function [T,P] = rev_isentropic(M,T0,P0)
gamma = 1.4;
T = T0./(1+M^2*(gamma-1)/2);
P = P0/(1+M^2*(gamma-1)/2)^(gamma/(gamma-1));
end

