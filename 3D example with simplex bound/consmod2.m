function [c, ceq] = consmod(input,fk,Xant,lambda,aux1,Hk,n0,b0,EU)


param

u1 = input(1);  u2 = input(2);
u3 = input(3);  

x0 = [u1 u2 u3]';

x1 = Xant(:,3);
x2 = Xant(:,2);
x3 = Xant(:,1);


% rho = 1/5*norm(x1-x2);

U = [x1-x0, x2-x0, x3-x0];
invU = inv(U);

r = 0.5*norm([(x1-x0)'*(x1-x0), (x2-x0)'*(x2-x0), (x3-x0)'*(x3-x0)]*invU);

Tr = L*r;

U2 = [x0 x1 x2 x3];

Ts = Etf(U2,L);

Lmin = lmin(U2);

Nl = 2*delt/Lmin;


% Constraints

c(1) = Ts + Nl - EU;
% c(1) = Ts + Nl - EU;

if aux1==1
    c(2) = b0 - n0'*x0;
else
    c(2) = n0'*x0 - b0;
end



ceq = [];