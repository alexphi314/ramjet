%philpott's graphing script

alt = [10 20 30 35];
T = [1345 288.9 60 26.27];
Isp = [6070 6106 6060 5991];

figure; hold on;
h1 = plot(T,alt,'DisplayName','Thrust (kN)');
h2 = plot(Isp,alt,'DisplayName','Isp (s-1)');
ylabel('Altitude (km)');
title('Altitude vs. Thrust and Isp');
obj = [h1 h2];
legend(obj,'Location','NortheastOutside');
print('altthrust','-djpeg');