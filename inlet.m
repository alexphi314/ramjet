function [P, T, P_0i, T_0i, u_i, M, comp_ratio, ds, eff, area] = inlet(M, P, T, m_dot, theta)
%Mach#, Ambient Pressure [Pa], Ambient Temperature [K], Desired Mass flow rate [kg/s], Row vector of deflections in degrees
%Ethen Daniels, ethend@umich.edu
%4/12/2017
%Calculates resulting flow conditions after multiple deflections

gamma = 1.4;
Cp = 1004.5;
R = 287;

%Initializing and conserving important quantities
P_inf = P;
T_inf = T;
M_inf = M;
rho_inf = P/(R*T);
rho = P/(R*T);
P_0inf = P*(1+(gamma-1)*M^2/2)^(gamma/(gamma-1));
T_0inf = T*(1+(gamma-1)*M^2/2);

%Calculates change accross a shock
for i = [1:size(theta,2)]   
    %Finds B for OS
    M_n=M;
    if theta(i) ~= 0
        fun = @(B) 2.*cot(B).*(M^2.*sin(B).^2 - 1)/(M^2.*(gamma+cos(B.*2))+2) - tan(theta(i)*pi/180);
        B=fsolve(fun,asin(1/M));
        M_n = M*sin(B);
    end

    %Calculates changes across shock
    P = P*(1 + 2*gamma/(gamma+1)*(M_n^2 -1));
    T = T*(1 + gamma*2/(gamma+1)*(M_n^2 -1))*(1 + (gamma-1)/2*M_n^2)/((gamma+1)/2*M_n^2);
    rho = rho*(gamma+1)*M_n^2/(2+(gamma-1)*M_n^2);
    M_n = sqrt((1+(gamma-1)/2*M_n^2)/(gamma*M_n^2-(gamma-1)/2));
    
    M = M_n;
    if theta(i)~=0
        M = M_n/sin(B-theta(i)*pi/180);
    end

end
%Assigns output values
u_i = M*sqrt(gamma*R*T);
P_0i = P*(1+(gamma-1)*M^2/2)^(gamma/(gamma-1));
T_0i = T*(1+(gamma-1)*M^2/2);
comp_ratio = rho/rho_inf;
ds = Cp*log(T/T_inf) - R*log(P/P_inf);
eff = P_0i/P_0inf;
area = m_dot/(rho*u_i);

end