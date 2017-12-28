function [M7,Astar,A7,T,Isp,L1,L2] = nozzle(P5,T5,M5,P1,mdot,mdotf,A4,w,u3,P3,A3,Pinlet,lengthinlet,thetaf,n)
    g = 1.4;
    R = 287;
    G = g.*(2./(g+1)).^((g+1)./2./(g-1));
    
    [T05, P05] = isentropic(M5,T5,P5);
    
    P7 = P1;
    P07 = P05;
    A5 = A4;
    
    M7 = ((2./(g-1)).*((P07./P7).^((g-1)./g)-1)).^0.5;
    a0 = (g.*R.*T05).^0.5;
    Astar = mdot*a0/(P05*G);
    
    A7 = Astar./M7.*(2./(g+1).*(1+(g-1)./2.*M7.^2)).^((g+1)./2./(g-1));
    
    
    [T7,~] = rev_isentropic(M7,T05,P05);
    a7 = (g.*R.*T7).^0.5;
    u7 = a7.*M7;
    
    %thrust calculation
    
    a = size(Pinlet);
    a = max(a);
    pai = 0;
    ai = 0;
    theta = [];
    for j = 1:n
        theta(j) = j.*thetaf./n;
    end
    %display(theta);
    for k = 2:a-1
            %display(theta(k-1));
            %display(lengthinlet(k-1));
            pa = Pinlet(k).*lengthinlet(k-1).*w.*sind(theta(k-1));
            pai = pai + pa;
            ai = lengthinlet(k-1).*w.*sind(theta(k-1)) + ai;
    end
    pai = pai + P3.*A3.*cosd(thetaf);
    ai = ai + A3.*cosd(thetaf);
    pae = P7.*A7;
    
    T = mdot.*(u3.*cosd(thetaf)-u7) + pai - pae - P1.*(ai - A7);
    display(T);
    Isp = T./mdotf./9.81;
    
    deltay = A5./w - Astar./w;
    L5t = deltay;
    
    deltay2 = A7./w - Astar./w;
    Lt7 = deltay2;
    
    L1 = L5t;
    L2 = Lt7;