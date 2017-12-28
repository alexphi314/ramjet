function [M5,T5,P5,pstar,p0star,Tstar] = flameholdera(M4,T4,P4,k)
    g = 1.4;
    R = 287;
    
    [T04,p04] = isentropic(M4,T4,P4);
    
    Tstar = T4.*((2 + (g-1).*M4.^2)./(g+1));
    pstar = P4.*M4.*(Tstar./T4).^0.5;
    p0star = p04.*M4.*((2 + (g-1).*M4.^2)./(g+1)).^((g+1)./2./(g-1));
    
    
    p05 = p04.*(1-k.*g.*M4.^(2)./2).*(1+(g-1)./2.*M4.^2).^(g./(1-g));
    
    syms M5;
    eqn1 = p05./p0star == 1./M5.*((g+1)./(2 + (g-1).*M5.^2)).^((g+1)./2./(g-1));
    soln = solve(eqn1,M5);
    a = size(soln);
    a = max(a);
    M5s = [];
    for k = 1:a
        M5s(k) = double(soln(k,1));
        t = isreal(M5s(k));
        if M5s(k) > 0 && t == 1
            if M4 < 1 && M5s(k) < 1
                M5 = M5s(k);
            elseif M4 > 1 && M5s(k) > 1
                M5 = M5s(k);
            end
        end
    end
    
    T05 = T04;
    
    [T5,P5] = rev_isentropic(M5, T05, p05);
    
end
    