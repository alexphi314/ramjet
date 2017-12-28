%inlet test script
%Generates the delta_theta for each section of the ramp
close all;
theta_final = 36; %in degrees
n = 5; %number of os
theta = [linspace(theta_final/n,theta_final/n,n),0]; %Determines the number ramps and the variation in theta
M = 3;
P = 1197;
T = 226;
m_dot = 100;
%Assigns values to vector based on number of steps

[P, T, P_0, T_0, u, M,B_vect, comp_ratio, ds, eff, area, s] = inlet_v(M, P, T, m_dot, theta,2);

%Gets overall dimensions of ramp
L=[];
H=[];
for  i = [1:size(s,2)]
    theta_t = sum(theta(1:i))*pi/180;
    L = [L, s(i)*cos(theta_t)];
    H = [H, s(i)*sin(theta_t)];
end
sum(L)
sum(H)