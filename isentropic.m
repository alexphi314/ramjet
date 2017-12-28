function [T0,P0] = isentropic(M,T,P)
gamma = 1.4;
T0 = T.*(1 + M.^2.*(gamma-1)./2);
P0 = P.*(1 + M.^2.*(gamma-1)./2).^(gamma./(gamma-1));
end

