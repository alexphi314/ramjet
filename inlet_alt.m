function [P, T, P_0, T_0, u, M,B_vect, comp_ratio, ds, eff, area, s,ramp_x,ramp_y] = inlet_alt(M, P, T, m_dot, theta, depth)
%Mach#, Ambient Pressure [Pa], Ambient Temperature [K], Desired Mass flow rate [kg/s], Row vector of deflections in degrees
%Ethen Daniels, ethend@umich.edu
%4/21/2017
%MODIFIED TO KEEP CONSTANT AREA
%Calculates resulting flow conditions after multiple deflections
%Output will be broken up into all involved regions. Quantities such as
%P, T, P_0, T_0, u, ds, and M will include data from each section of the flow,
%from ambient to inlet, including regions between oblique shocks. B_vect is
%the angle Beta of the OS relative to the normal flow in that section.
%The elements are ordered in the sequence the occur moving from ambient to
%inlet. All output units are standard SI.

gamma = 1.4;
Cp = 1004.5;
R = 287;

%Initializes vectors
P = [P,zeros(1,size(theta,2))];
T = [T,zeros(1,size(theta,2))];
M = [M,zeros(1,size(theta,2))];
rho = [P(1)/(R*T(1)),zeros(1,size(theta,2))];
P_0 = [P(1)*(1+(gamma-1)*M(1)^2/2)^(gamma/(gamma-1)),zeros(1,size(theta,2))];
T_0 = [T(1)*(1+(gamma-1)*M(1)^2/2),zeros(1,size(theta,2))];
ds = zeros(1,size(theta,2)+1);
B_vect = zeros(1,size(theta,2)-1);

%Calculates the changes in each section of the flow
for i = [1:size(theta,2)]   
    %Finds B for OS
    M_n=M(i);
    if theta(i) ~= 0
        fun = @(B) 2.*cot(B).*(M(i)^2.*sin(B).^2 - 1)/(M(i)^2.*(gamma+cos(B.*2))+2) - tan(theta(i)*pi/180);
        B=fsolve(fun,asin(1/M(i)));
        M_n = M(i)*sin(B);
    end

    %Calculates changes across shock
    P(i+1) = P(i)*(1 + 2*gamma/(gamma+1)*(M_n^2 -1));
    T(i+1) = T(i)*(1 + gamma*2/(gamma+1)*(M_n^2 -1))*(1 + (gamma-1)/2*M_n^2)/((gamma+1)/2*M_n^2);
    rho(i+1) = P(i+1)/(R*T(i+1));
    M_n = sqrt((1+(gamma-1)/2*M_n^2)/(gamma*M_n^2-(gamma-1)/2));
    
    %Assigns new Mach number bases on whether there was an OS or NS
    M(i+1) = M_n;
    if theta(i)~=0
        M(i+1) = M_n/sin(B-theta(i)*pi/180);
    end
    
    %Calculates the change in entropy for current shock
    ds(i+1) = Cp*log(T(i+1)/T(i)) - R*log(P(i+1)/P(i));
    
    %Saves B in a vector for output
    B_vect(i)=B;
end
B_vect = B_vect(1:end-1);
%Assigns output values
u = M.*sqrt(gamma*R.*T);
P_0 = P.*(1+(gamma-1).*M.^2/2).^(gamma/(gamma-1));
T_0 = T.*(1+(gamma-1).*M.^2/2);
comp_ratio = rho(end)/rho(1);
eff = P_0(end)/P_0(1);
area = 1.7814;
[s,ramp_x,ramp_y] = ramp(B_vect, theta, area, depth);
end