%%%%%%%%%%%%%%%%%%%%% Computation of Ea %%%%%%%%%%%%%%%%%%%%%
% Computes the norm of the truncation gradient error of function f
% evaluated a vopt.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ea = Etf(Udata,D2)

[nu,b] = size(Udata);

if not(nu == b-1)
    disp('Error: Matrix dimensions are incorrect')
    return
end

switch nu

   case 2

       x1s = Udata(:,1);  % u_{k-1}
       x2s = Udata(:,2);  % u_k
       x3s = Udata(:,3);  % u

       U = [x3s-x2s, x3s-x1s];
       q = [(x3s-x2s)'*(x3s-x2s), (x3s-x1s)'*(x3s-x1s)];
%        Et(i,j) = dd2/2*norm(q*inv(U));
        
       w = -D2/2*q*inv(U);
                
       H = [x3s x2s x1s]'*[x3s x2s x1s];
       f = [-(1/D2*w + x3s')*[x3s x2s x1s]]';
       Aeq = [1 1 1];
       beq = 1;
       lb = [0 0 0]';
        
       d0 = [0.3 0.3 0.4];
        
        
       optionsx = optimset('Display','off','LargeScale','off','MaxFunEvals',500,'MaxIter',300,...
                           'TolFun',1e-10,'TolX',1e-10,'TolCon',1e-9);
    
       [ds,fval] = quadprog(H,f,[],[],Aeq,beq,lb,[],d0,optionsx);
        
       dopt = [ds(1) ds(2) ds(3)]';
        
       vopt = [x3s x2s x1s]*dopt;
        
       efopt = -D2/2*q*inv(U) + D2*(x3s-vopt)';
        
       Ea = norm(efopt);

   case 3

       x1s = Udata(:,1);  % u_{k-2}
       x2s = Udata(:,2);  % u_{k-1}
       x3s = Udata(:,3);  % u_k
       x4s = Udata(:,4);  % u

       U = [x4s-x3s, x4s-x2s, x4s-x1s];
       q = [(x4s-x3s)'*(x4s-x3s), (x4s-x2s)'*(x4s-x2s), (x4s-x1s)'*(x4s-x1s)];
%        Et(i,j) = dd2/2*norm(q*inv(U));
        
       w = -D2/2*q*inv(U);
                
       H = [x4s x3s x2s x1s]'*[x4s x3s x2s x1s];
       f = [-(1/D2*w + x4s')*[x4s x3s x2s x1s]]';
       Aeq = [1 1 1 1];
       beq = 1;
       lb = [0 0 0 0]';
        
       d0 = [0.25 0.25 0.25 0.25];
        
        
       optionsx = optimset('Display','off','LargeScale','off','MaxFunEvals',500,'MaxIter',300,...
                           'TolFun',1e-10,'TolX',1e-10,'TolCon',1e-9);
    
       [ds,fval] = quadprog(H,f,[],[],Aeq,beq,lb,[],d0,optionsx);
        
       dopt = [ds(1) ds(2) ds(3) ds(4)]';
        
       vopt = [x4s x3s x2s x1s]*dopt;
        
       efopt = -D2/2*q*inv(U) + D2*(x4s-vopt)';
        
       Ea = norm(efopt);

    otherwise

      disp('Only 2 and 3-dimensional input spaces are supported')

end

