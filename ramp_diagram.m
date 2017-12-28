function [] = ramp_diagram(s,theta,B_vect)
%Generates Ramp geometry along with OS geometry
%Ethen Daniels, 4/19/2017, ethend@umich.edu
close ALL;

x = zeros(1,size(theta,2));
y = zeros(1,size(theta,2));
for i = [2:size(theta,2)]
    x(i) = x(i-1) + s(i-1)*cos(sum(theta(1:i-1))*pi/180);
    y(i) = y(i-1) + s(i-1)*sin(sum(theta(1:i-1))*pi/180);
end
figure(1);
hold;
plot(x,y,'-k','LineWidth',1);
axis([0 10 0 10]);
end