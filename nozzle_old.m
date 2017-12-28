function [M7,Astar,A7,T,Isp,L1,L2] = nozzle(P5,T5,M5,P1,mdot,mdotf,A4,w)
    g = 1.4;
    R = 287;
    G = g.*(2./(g+1)).^((g+1)./2./(g-1));
    
    [T05, P05] = isentropic(M5,T5,P5);
    
    P7 = P1;
    P07 = P05;
    A5 = A4;
    
    %M7 = ((2./(g-1)).*((P07./P7).^((g-1)./g)-1)).^0.5;
    a0 = (g.*R.*T05).^0.5;
    Astar = mdot.*a0./P07./G;
    A7 = A5;
    
    syms M7;
    eqn1 = A7 == Astar./M7.*(2./(g+1).*(1+(g-1)./2.*M7.^2)).^((g+1)./2./(g-1));
    soln = solve(eqn1,M7);
    a = size(soln);
    a = max(a);
    for k = 1:a
    M7s(k) = double(soln(k,1));
    t = isreal(M7s(k));
    if M7s(k) > 0 && t==1
        if M7s(k) > 1
            M7 = M7s(k);
%         elseif M3 >1 && M4s(k) >1
%             M4 = M4s(k);
        end
    end
end
    
    
    [T7,P7] = rev_isentropic(M7,T05,P05);
    a7 = (g.*R.*T7).^0.5;
    u7 = a7.*M7;
    %display(u7);
    
    T = mdot.*u7 + A7.*(P7-P1);
    %display(P7-P1);
    Isp = T./mdotf./g;
    
    deltay = A5./w - Astar./w;
    L5t = deltay;
    
    deltay2 = A7./w - Astar./w;
    Lt7 = deltay2;
    
    L1 = L5t;
    L2 = Lt7;
end