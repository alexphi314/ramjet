function [M4,T4,P4,Astar] = diffusera(M3,T3,P3,A4,mdot)
g = 1.4;
R = 287;
Cp = 1004.5;
G = g.*(2./(g+1)).^((g+1)./2./(g-1));
%display(G);

[T03,p03] = isentropic(M3,T3,P3);
a0 = (g.*R.*T03).^0.5;

Astar = 1.6149;
display(Astar);

syms M4;
eqn1 = A4./Astar == 1./M4.*(2./(g+1).*(1+(g-1)./2.*M4.^2)).^((g+1)./2./(g-1));
soln = solve(eqn1,M4);

a = size(soln);
a = max(a);
M4s = [];
for k = 1:a
    M4s(k) = double(soln(k,1));
    t = isreal(M4s(k));
    if M4s(k) > 0 && t==1
        if M4s(k) < 1
            M4 = M4s(k);
%         elseif M3 >1 && M4s(k) >1
%             M4 = M4s(k);
        end
    end
end

[T4,P4] = rev_isentropic(M4,T03,p03);

end






