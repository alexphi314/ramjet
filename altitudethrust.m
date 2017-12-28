%find resulting thrust at different altitudes

%% Alex Philpott Ethen Daniels Dylan Ma
clear;
axis equal;
clc;
close all;
R = 287;
M1 = 3;
qhv = 120e6;
g = 1.4;
mdot = 100;
T5 = 1800;
Cp = 1004.5;
Cv = 717.5;
w = 2;

H = 30;
if H == 30
    P1 = 1197;
    T1 = 226;
    rho1 = P1./R./T1;
elseif H == 20
    P1 = 5475.21;
    T1 = 216.649;
    rho1 = 0.08803;
elseif H == 10
    P1 = 26437.3;
    T1 = 223.150;
    rho1 = 0.41271;
elseif H == 35
    P1 = 558.969;
    T1 = 237.048;
    rho1 = 0.00821;
end
    

[T01,P01] = isentropic(M1,T1,P1);

Mplot = [];
Pplot = [];
P0plot = [];
Tplot = [];
T0plot = [];
uplot = [];
hplot = [];
splot = [];
Lplot = [];
Thrs = [];
Ltrack = [];
metals = [];

%for k = 1:100
A4 = 6; %+ k./10;


%% Inlet
thetaf = 36;
n = 5; %number of oblique shocks
theta = [];
for k = 1:n
    theta(1,k) = thetaf./n;
end
theta(end+1) = 0;

[P,T,P_0,T_0,u,M,B_vect,prat,ds,eff,A3,length,ramp_x,ramp_y] = inlet_v(M1,P1,T1,mdot,theta,w);
Pinlet = P;
lengthinlet = length;
P3 = P(end);
T3 = T(end);
P03 = P_0(end);
T03 = T_0(end);
u3 = u(end);
M3 = M(end);
a3 = (g.*R.*T3).^0.5;
c = size(M);
Ms = [];
Ps = [];
P0s = [];
Ts = [];
T0s = [];
us = [];
ss = [];
hs = [];
Ls = [];
Linlet = 0;

for k = 1:max(c)
    Ms(2.*k-1) = M(k);
    Ms(2.*k) = M(k);
    Ps(2.*k-1) = P(k);
    Ps(2.*k) = P(k);
    P0s(2.*k-1) = P_0(k);
    P0s(2.*k) = P_0(k);
    Ts(2.*k-1) = T(k);
    Ts(2.*k) = T(k);
    T0s(2.*k-1) = T_0(k);
    T0s(2.*k) = T_0(k);
    us(2.*k-1) = u(k);
    us(2.*k) = u(k);
    ss(2.*k-1) = ds(k);
    ss(2.*k) = ds(k);
    h = Cp.*T(k);
    hs(2.*k-1) = h;
    hs(2.*k) = h;
end
for k = 1:n
    L = length(1,k).*cosd(k.*theta(1,1));
    Ls(2.*k) = Linlet;
    Ls(2.*k+1) = Linlet + L;
    Linlet = Linlet + L;
end
Ls = Ls(2:end);
Mplot = [Mplot Ms];
Pplot = [Pplot Ps];
Tplot = [Tplot Ts];
P0plot = [P0plot P0s];
T0plot = [T0plot T0s];
uplot = [uplot us];
splot = [splot ss];
hplot = [hplot hs];
Lplot = [Lplot Ls];
Mplot = Mplot(2:end-1);
Pplot = Pplot(2:end-1);
Tplot = Tplot(2:end-1);
P0plot = P0plot(2:end-1);
T0plot = T0plot(2:end-1);
uplot = uplot(2:end-1);
splot = splot(2:end-1);
hplot = hplot(2:end-1);
Lplot(end+1) = Linlet;
Lplot(end+1) = Linlet;

%% Diffuser

L23 = 0.5.*cosd(thetaf);
Lplot(end) = Lplot(end) + L23;
[M4,T4,P4,Astar] = diffuser(M3,T3,P3,A4,mdot);
rho4 = P4./R./T4;
a4 = (g.*R.*T4).^0.5;
u4 = M4.*a4;
Ldiff = A4./w - A3./w;

[T04, P04] = isentropic(M4,T4,P4);

Ms = linspace(M3,M4);
A = Astar./Ms.*(2./(g+1).*(1+(g-1)./2.*Ms.^2)).^((g+1)./2./(g-1));
p0 = P03.*ones(1,100);
T0 = T03.*ones(1,100);
T = T0./(1+(g-1)./2.*Ms.^2);
P = p0./(1+(g-1)./2.*Ms.^2).^(g./(g-1));
a = (g.*R.*T).^0.5;
u = a.*Ms;
Ls = linspace(Lplot(end),Lplot(end)+Ldiff);

figure; hold on;
plot(Ls,-A./w,'b');
diff_x = Ls;
diff_y = -A./w;
plot([Ls(1) Ls(1)],[-A(1)./w 0],'b');
plot([Ls(1) Ls(end)],[0 0],'b');
plot([Ls(end) Ls(end)],[-A(end)./w 0],'b');
xlabel('Position (m)');
ylabel('Height (m)');
ylim([-3.5,0.5]);
xlim([7 10]);
title('Diagram of Diffuser');
print('diffuser','-djpeg');

ds = zeros(1,100);
ds = ds + splot(end);
h1 = Cp.*T3;
u1 = a3.*M3;
h = h1 + 0.5.*u1.^2 - 0.5.*u.^2;

Mplot = [Mplot Ms];
Pplot = [Pplot P];
Tplot = [Tplot T];
P0plot = [P0plot p0];
T0plot = [T0plot T0];
uplot = [uplot u];
splot = [splot ds];
hplot = [hplot h];
Lplot = [Lplot Ls];

%% Combustor
[M5,M41,T41,P41,P5,T5,rho41,rho5,f,Lcc,pstarf,p0starf,Tstarf,Pstarc,Tstarc] = combustionsection(M4,T4,P4,k,qhv,T5,A4);
mdotf = mdot.*f;
a41 = (g.*R.*T41).^0.5;
u41 = a41.*M41;
[T0starc,P0starc] = isentropic(1,Tstarc,Pstarc);
[T041, P041] = isentropic(M41,T41,P41);

%flameholder
G = rho4.*u4;
h04 = Cp.*T04;
rhos = linspace(rho4,rho41);

M = (2.*pstarf.^2./((g+1).*(rhos.^2.*R.^2.*Tstarf.^2)-(g-1).*pstarf.^2)).^0.5;
p = pstarf./M.*((g+1)./(2+(g-1).*M.^2)).^0.5;
p0 = p0starf./M.*((g+1)./(2+(g-1).*M.^2)).^(-(g+1)./2./(g-1));
T = Tstarf.*((g+1)./(2+(g-1).*M.^2));
T0 = T04.*ones(1,100);
a = (g.*R.*T).^0.5;
u = M.*a;

a5 = (g.*R.*T5).^0.5;
u5 = M5.*a5;

h = h04 - G.^2./2./rhos.^2;
s = Cv.*log(p./P4)-Cp.*log(rhos./rho4);
s = s + splot(end);
Lplot(end) = Lplot(end) + 0.4;
L = linspace(Lplot(end),Lplot(end)+0.4);

Mplot = [Mplot M];
Pplot = [Pplot p];
Tplot = [Tplot T];
P0plot = [P0plot p0];
T0plot = [T0plot T0];
uplot = [uplot u];
splot = [splot s];
hplot = [hplot h];
Lplot = [Lplot L];

%combustion chamber
rhos = linspace(rho41,rho5);
F = P41 + rho41.*u41.^2;
G = rho41.*u41;
h = g./(g-1)./rhos.*(F - G.^2./rhos);
s = -Cp.*log(rhos./rho41) + Cv.*log(F./P41 - G.^2./P41./rhos);
s = s + splot(end);

Ms = (Pstarc./(rhos.*R.*Tstarc.*(g+1)-g.*Pstarc)).^0.5;
p = Pstarc.*(1+g)./(1+g.*Ms.^2);
T = Tstarc.*Ms.^2.*(1+g).^2./(1+g.*Ms.^2).^2;
p0 = P0starc.*p./Pstarc.*((2 + (g-1).*Ms.^2)./(g+1)).^(g./(g-1));
T0 = T.*T0starc./Tstarc.*((2 + (g-1).*Ms.^2)./(g+1));
a = (g.*R.*T).^0.5;
u = a.*Ms;

Mplot = [Mplot Ms];
Pplot = [Pplot p];
Tplot = [Tplot T];
P0plot = [P0plot p0];
T0plot = [T0plot T0];
uplot = [uplot u];
splot = [splot s];
hplot = [hplot h];

Lcomb = Lcc;
L = linspace(Lplot(end),Lplot(end)+Lcomb);
comb_x = L;
comb_y = -A4./w;
Lplot = [Lplot L];
%plot([Lplot(end)-Lcomb-0.2,Lplot(end)],[-A4./w -A4./w],'b');
%plot([Lplot(end)-Lcomb-0.2,Lplot(end)],[0 0],'b');

%% Nozzle

[M7,Astar,A7,Thr,Isp,Lnoz1,Lnoz2] = nozzle(P5,T5,M5,P1,mdot,mdotf,A4,w,u3,P3,A3,Pinlet,lengthinlet,thetaf,n);
[T05,P05] = isentropic(M5,T5,P5);
A5 = A4;

Ms = linspace(M5,1);
Ms = [Ms linspace(1,M7)];
p0 = P05.*ones(1,200);
T0 = T05.*ones(1,200);
T = T0./(1+(g-1)./2.*Ms.^2);
P = p0./(1+(g-1)./2.*Ms.^2).^(g./(g-1));
a = (g.*R.*T).^0.5;
u = a.*Ms;
A = (A5./Ms.*(2./(g+1).*(1+(g-1)./2.*Ms.^2)).^((g+1)./2./(g-1)))./(1./M5.*(2./(g+1).*(1+(g-1)./2.*M5.^2)).^((g+1)./2./(g-1)));
L = linspace(Lplot(end),Lplot(end)+Lnoz1);
Lthroat = L(end);
L = [L linspace(L(end),L(end) + Lnoz2)];

P6 = P(100);
T6 = T(100);
T06 = T0(100);
P06 = p0(100);
u6 = u(100);
Athroat = A(100);
M6 = 1;
rho6 = P6./T6./R;

figure; hold on;
plot(L,-A./w,'b');
plot([L(1) L(1)],[-A(1)./w 0],'b');
%plot([L(1) L(end)],[0 0],'b');

ds = zeros(1,200);
ds = ds + splot(end);
h1 = Cp.*T5;
a5 = (g.*R.*T5).^0.5;
u1 = a5.*M5;
h = h1 + 0.5.*u1.^2 - 0.5.*u.^2;

Mplot = [Mplot Ms];
Pplot = [Pplot P];
Tplot = [Tplot T];
P0plot = [P0plot p0];
T0plot = [T0plot T0];
uplot = [uplot u];
splot = [splot ds];
hplot = [hplot h];
Lplot = [Lplot L];

% Ms = linspace(1,M7);
% p0 = P05.*ones(1,100);
% T0 = T05.*ones(1,100);
% T = T0./(1+(g-1)./2.*Ms.^2);
% P = p0./(1+(g-1)./2.*Ms.^2).^(g./(g-1));
% a = (g.*R.*T).^0.5;
% u = a.*Ms;
% A = Astar./Ms.*(2./(g+1).*(1+(g-1)./2.*Ms.^2)).^((g+1)./2./(g-1));
% L = linspace(Lplot(end),Lplot(end)+Lnoz2);

M7 = Ms(end);
P07 = p0(end);
T07 = T0(end);
T7 = T(end);
P7 = P(end);
u7 = u(end);
rho7 = P7./T7./R;

% plot(L,-A./w,'b');
% noz2_x = L;
% noz2_y = -A./w;
plot([L(end) L(end)],[-A(end)./w 0],'b');
plot([L(1) L(end)],[0 0],'b');
xlabel('Position (m)');
ylabel('Height (m)');
title('Nozzle Diagram');
ylim([-6 0.5]);
xlim([26.5 33]);
print('nozzle','-djpeg');

ds = zeros(1,100);
ds = ds + splot(end);
h = h1 + 0.5.*u1.^2 - 0.5.*u.^2;

% Mplot = [Mplot Ms];
% Pplot = [Pplot P];
% Tplot = [Tplot T];
% P0plot = [P0plot p0];
% T0plot = [T0plot T0];
% uplot = [uplot u];
% splot = [splot ds];
% hplot = [hplot h];
% Lplot = [Lplot L];

L = Lplot(end);
metal = 2.*(w + A4./w).*Lcomb;
metals = [metals metal];
Ltrack = [Ltrack L];
Thrs = [Thrs Thr];
%end

%% Plots

%All properties across engine

figure; hold on;
h1 = plot(Lplot,Mplot,'DisplayName','Mach Number');
h2 = plot(Lplot,Pplot./P1,'DisplayName','Pressure');
h3 = plot(Lplot,P0plot./P01,'DisplayName','Stagnation Pressure');
h4 = plot(Lplot,Tplot./T1,'DisplayName','Temperature');
h5 = plot(Lplot,T0plot./T01,'DisplayName','Stagnation Temperature');
h6 = plot(Lplot,uplot./u1,'DisplayName','Speed');
xlabel('Length (m)');
ylabel('Normalized Properties');
title('Properties vs. Length');
obj = [h1 h2 h3 h4 h5 h6];
legend(obj,'Location','NortheastOutside');
print('properties','-djpeg');

figure; hold off;
plot(splot,hplot);
xlabel('Change in Entropy (J/kg-K)');
ylabel('Enthalpy (J/kg)');
title('Mollier Diagram');
print('hs','-djpeg');

figure;
plot(Lplot,hplot);