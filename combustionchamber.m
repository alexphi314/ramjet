function [M6,P6,T6,f,pstar,Tstar] = combustionchamber(M5,T5,P5,qhv,Tcomb)
    g = 1.4;
    R = 287;
    Cp = 1004.5;
    
    [T05,~] = isentropic(M5,T5,P5);
    
    Tstar = T5./M5.^2.*((1 + g.*M5.^2)./(g + 1)).^2;
    display(Tstar);
    pstar = P5.*(1 + g.*M5.^2)./(1 + g);
    
    if Tcomb > Tstar
%         M6 = 1;
%         P6 = pstar.*(1+g)./(1+g.*M6.^2);
%         T6 = Tstar;
    syms M6;
    T6 = 0.95.*Tstar;
    eqn1 = T6./Tstar == M6.^2.*((1 + g)./(1 + g.*M6.^2)).^2;
    soln = solve(eqn1,M6);
    a = size(soln);
    a = max(a);
    M6s = [];
    for k = 1:a
        M6s(k) = double(soln(k,1));
        t = isreal(M6s(k));
        if M6s(k) > 0 && t==1
            if M5 < 1 && M6s(k) < 1
                M6 = M6s(k);
%             elseif M5 >1 && M6s(k) >1
%                 M6 = M6s(k);
            end
        end
    end
    
    P6 = pstar.*(1+g)./(1+g.*M6.^2);
    elseif Tcomb < Tstar
    syms M6;
    eqn1 = Tcomb./Tstar == M6.^2.*((1 + g)./(1 + g.*M6.^2)).^2;
    soln = solve(eqn1,M6);
    a = size(soln);
    a = max(a);
    M6s = [];
    for k = 1:a
        M6s(k) = double(soln(k,1));
        t = isreal(M6s(k));
        if M6s(k) > 0 && t==1
            if M5 < 1 && M6s(k) < 1
                M6 = M6s(k);
            elseif M5 >1 && M6s(k) >1
                M6 = M6s(k);
            end
        end
    end
    
    P6 = pstar.*(1+g)./(1+g.*M6.^2);
    T6 = Tcomb;
    end
    
    [T06,p06] = isentropic(M6,T6,P6);
    
    q = Cp.*(T06 - T05);
    f = q./qhv;
end