function [c, ceq] = consmod(input,fk,Xant,lambda,aux1,Hk,n0,b0)


param

u1 = input(1);
u2 = input(2);

x0 = [u1 u2]';

x1 = Xant(:,2);
x2 = Xant(:,1);

% rho = 1/5*norm(x1-x2);

% U = [x1-x0, x2-x0];
% invU = inv(U);
% 
% r = 0.5*norm([(x1-x0)'*(x1-x0), (x2-x0)'*(x2-x0)]*invU);
% 
% Tr = L*r;

U2 = [x0 x1 x2];

Ts = Etf(U2,L);

Lmin = lmin(U2);

Nl = 2*delt/Lmin;


% Constraints

% c(1) = Ts - EU;
c(1) = Ts + Nl - EU;

if aux1==1
    c(2) = b0 - n0'*x2;
else
    c(2) = n0'*x2 - b0;
end



ceq = [];