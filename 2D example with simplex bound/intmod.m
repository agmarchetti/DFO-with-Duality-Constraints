
function [f, g] = intmod(input,fk,Xant,lambda,aux1,Hk,n0,b0)


u1 = input(1);
u2 = input(2);

u = [u1 u2]';
uk = Xant(:,1);

q = 0.5*(u-uk)'*Hk*(u-uk);

m = q + fk + lambda*(u - uk);

f = m;

g = [];