%Ethen Daniels, ethend@umich.edu
%4/12/2017
%Calculates inlet properties for multiple geometries using inlet.m
close ALL;
clear;
%Inputs
P = 1197;   %Pressure at 30km
T = 226;    %Temperature at 30km
M = 3;
m_dot = 100;
geo = [0:10]; %Number of increments till desired theta;


%Initializes figures and vectors
for i = [1:10]
    figure(i);
    hold;
end
for theta_final = [35 45]
P_i = zeros(1,size(geo,2));
T_i = zeros(1,size(geo,2));
P_0i = zeros(1,size(geo,2));
T_0i = zeros(1,size(geo,2));
u_i = zeros(1,size(geo,2));
M_i = zeros(1,size(geo,2));
comp_ratio = zeros(1,size(geo,2));
ds = zeros(1,size(geo,2));
eff = zeros(1,size(geo,2));
area = zeros(1,size(geo,2));
for n = geo
    %Generates the delta_theta for each section of the ramp
    theta = n;
    if n>0
        theta = [linspace(theta_final/n,theta_final/n,n),0];
    end    
    
    %Assigns values to vector based on number of steps
    [P_i(n+1), T_i(n+1), P_0i(n+1), T_0i(n+1), u_i(n+1), M_i(n+1), comp_ratio(n+1), ds(n+1), eff(n+1), area(n+1)] = inlet(M, P, T, m_dot, theta);

end

figure(1);
plot(geo,P_i,'x');
title('Inlet Pressure');
xlabel('Number of steps');
ylabel('Pressure [Pa]');

figure(2);
plot(geo,T_i,'x');
title('Inlet Temperature');
xlabel('Number of steps');
ylabel('Temperature [K]');

figure(3);
plot(geo,P_0i,'x');
title('Inlet Stagnation Pressure');
xlabel('Number of steps');
ylabel('Pressure [Pa]');

figure(4);
plot(geo,T_0i,'x');
title('Inlet Stagnation Temperature');
xlabel('Number of steps');
ylabel('Temperature [K]');

figure(5);
plot(geo,u_i,'x');
title('Inlet Flow Velocity');
xlabel('Number of steps');
ylabel('Velocity [m/s]');


figure(6);
plot(geo,M_i,'x');
title('Inlet Flow Mach #');
xlabel('Number of steps');
ylabel('Mach #');


figure(7);
plot(geo,comp_ratio,'x');
title('Inlet Compression Ratio');
xlabel('Number of steps');
ylabel('Compression [rho_i/rho_i_n_f]');

figure(8);
plot(geo,ds,'x');
title('Ds accross Inlet');
xlabel('Number of steps');
ylabel('Ds [J/(kgK)]');

figure(9);
plot(geo,eff,'x');
title('Inlet Efficiency');
xlabel('Number of steps');
ylabel('Stagnation Pressure Ratio [P_0_i/P_0_i_n_f]');

figure(10);
plot(geo,area,'x');
title('Inlet Area');
xlabel('Number of steps');
ylabel('Area [m^2]');

end
for i = [1:10]
    figure(i);
    legend('25 deg','30 deg','35 deg','Location','NorthEastOutside');
end