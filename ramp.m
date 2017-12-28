function [s,ramp_x,ramp_y] = ramp(B_vect, theta, area, depth)
%Generates Ramp geometry based on data from inlet_v.m
%Causes OS to converge to leading endge of cowl
%Ethen Daniels, 4/18/17, ethend@umich.edu

%See Attatched diagram and derivation to understand notation
theta = theta*pi/180;
s = zeros(1,size(B_vect,2));
c = zeros(1,size(B_vect,2));
alpha = B_vect - theta(1);
H = area/depth; %Length of NS at inlet
s(end) = H/tan(alpha(end));
c(end) = H/sin(alpha(end));

%Solves for ramp and oblique shock geometry using method derived (See
%attatched)
for i = linspace(size(B_vect,2)-1,1,size(B_vect,2)-1)
    c(i) = c(i+1)*sin(pi - alpha(i+1)-theta(1))/sin(alpha(i));
    s(i) = c(i)*cos(alpha(i)) - sqrt(c(i+1)^2 - c(i)^2*sin(alpha(i))^2);
end

%Generates plot data for geometry
ramp_x = zeros(1,size(theta,2));
ramp_y = zeros(1,size(theta,2));
OS_x = zeros(size(theta,2)-1,2);
OS_y = zeros(size(theta,2)-1,2);
for i = [2:size(theta,2)]
    ramp_x(i) = ramp_x(i-1) + s(i-1)*cos(sum(theta(1:i-1)));
    ramp_y(i) = ramp_y(i-1) + s(i-1)*sin(sum(theta(1:i-1)));
    OS_x(i-1,1) = ramp_x(i-1);
    OS_x(i-1,2) = OS_x(i-1,1) + c(i-1)*cos(alpha(i-1)+sum(theta(1:i-1)));
    OS_y(i-1,1) = ramp_y(i-1);
    OS_y(i-1,2) = OS_y(i-1,1) + c(i-1)*sin(alpha(i-1)+sum(theta(1:i-1)));
end
figure;
hold;
for i = [1:size(OS_x,1)]
plot(OS_x(i,:),OS_y(i,:),'-b');
end
plot(ramp_x,ramp_y,'-k','LineWidth',1);
title('Inlet Ramp and OS Geometry');
xlabel('Length [m]');
ylabel('Height [m]');
legend('Oblique Shock','Location','NorthEast');
axis([0 ceil(max(ramp_x)) 0 ceil(max(ramp_x))]);

end
