
function [f, g] = intmod(input,fk,Xant,lambda,aux1,Hk,n0,b0,EU)


u1 = input(1);
u2 = input(2);
u3 = input(3);

u = [u1 u2 u3]';
uk = Xant(:,3);

q = 0.5*(u-uk)'*Hk*(u-uk);

m = q + fk + lambda*(u - uk);

f = m;

g = [];