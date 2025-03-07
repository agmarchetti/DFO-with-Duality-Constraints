%%%%%%%%%%%%%%%%%%%%% Computation of lmin %%%%%%%%%%%%%%%%%%%%%
% Computes the minimum distance between complement affine subspaces 
% generated from the nu+1 points in matrix U.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Lmin = lmin(U)

[nu,b] = size(U);

if not(nu == b-1)
    disp('Error: Matrix dimensions are incorrect')
    return
end

switch nu

   case 2

      % 3 point to line distances 

      n1 = null(U(:,2)' - U(:,3)');
      l1 = abs(n1'*(U(:,1) - U(:,2)))/norm(n1);

      n2 = null(U(:,1)' - U(:,3)');
      l2 = abs(n2'*(U(:,2) - U(:,1)))/norm(n2);

      n3 = null(U(:,1)' - U(:,2)');
      l3 = abs(n3'*(U(:,3) - U(:,1)))/norm(n3);

      Lmin = min([l1,l2,l3]);

   case 3

      % 4 point to plane distances

      R = [U(:,1) - U(:,2), U(:,1) - U(:,3)];
      n1 = null(R');
      l1 = abs(n1'*(U(:,4) - U(:,1)))/norm(n1);

      R = [U(:,1) - U(:,2), U(:,1) - U(:,4)];
      n2 = null(R');
      l2 = abs(n2'*(U(:,3) - U(:,1)))/norm(n2);

      R = [U(:,1) - U(:,3), U(:,1) - U(:,4)];
      n3 = null(R');
      l3 = abs(n3'*(U(:,2) - U(:,1)))/norm(n3);

      R = [U(:,2) - U(:,3), U(:,2) - U(:,4)];
      n4 = null(R');
      l4 = abs(n4'*(U(:,1) - U(:,2)))/norm(n4);

      % 3 line to line distances

      R = [U(:,1) - U(:,2), U(:,3) - U(:,4)];
      n5 = null(R');
      l5 = abs(n5'*(U(:,1) - U(:,3)))/norm(n5);

      R = [U(:,1) - U(:,3), U(:,2) - U(:,4)];
      n6 = null(R');
      l6 = abs(n6'*(U(:,1) - U(:,2)))/norm(n6);

      R = [U(:,2) - U(:,3), U(:,1) - U(:,4)];
      n7 = null(R');
      l7 = abs(n7'*(U(:,2) - U(:,4)))/norm(n7);

      % Computation of lmin
      
%       [l1 l2 l3 l4 l5 l6 l7]

      Lmin = min([l1 l2 l3 l4 l5 l6 l7]);

    case 4

      % 5 point to hyperplane distances

      R = [U(:,2) - U(:,3), U(:,2) - U(:,4), U(:,2) - U(:,5)];
      n1 = null(R');
      l1 = abs(n1'*(U(:,1) - U(:,2)))/norm(n1);

      R = [U(:,1) - U(:,3), U(:,1) - U(:,4), U(:,1) - U(:,5)];
      n2 = null(R');
      l2 = abs(n2'*(U(:,2) - U(:,3)))/norm(n2);

      R = [U(:,1) - U(:,2), U(:,1) - U(:,4), U(:,1) - U(:,5)];
      n3 = null(R');
      l3 = abs(n3'*(U(:,3) - U(:,1)))/norm(n3);

      R = [U(:,1) - U(:,2), U(:,1) - U(:,3), U(:,1) - U(:,5)];
      n4 = null(R');
      l4 = abs(n4'*(U(:,4) - U(:,1)))/norm(n4);

      R = [U(:,1) - U(:,2), U(:,1) - U(:,3), U(:,1) - U(:,4)];
      n5 = null(R');
      l5 = abs(n5'*(U(:,5) - U(:,1)))/norm(n5);

      % 10 line to plane distances

      R = [U(:,1) - U(:,2), U(:,3) - U(:,4), U(:,3) - U(:,5)];
      n6 = null(R');
      l6 = abs(n6'*(U(:,1) - U(:,3)))/norm(n6);

      R = [U(:,2) - U(:,3), U(:,1) - U(:,4), U(:,1) - U(:,5)];
      n7 = null(R');
      l7 = abs(n7'*(U(:,2) - U(:,1)))/norm(n7);

      R = [U(:,3) - U(:,4), U(:,1) - U(:,2), U(:,1) - U(:,5)];
      n8 = null(R');
      l8 = abs(n8'*(U(:,3) - U(:,1)))/norm(n8);

      R = [U(:,4) - U(:,5), U(:,1) - U(:,2), U(:,1) - U(:,3)];
      n9 = null(R');
      l9 = abs(n9'*(U(:,4) - U(:,1)))/norm(n9);

      R = [U(:,1) - U(:,3), U(:,2) - U(:,4), U(:,2) - U(:,5)];
      n10 = null(R');
      l10 = abs(n10'*(U(:,1) - U(:,2)))/norm(n10);

      R = [U(:,1) - U(:,4), U(:,2) - U(:,3), U(:,2) - U(:,5)];
      n11 = null(R');
      l11 = abs(n11'*(U(:,1) - U(:,2)))/norm(n11);

      R = [U(:,1) - U(:,5), U(:,2) - U(:,3), U(:,2) - U(:,4)];
      n12 = null(R');
      l12 = abs(n12'*(U(:,1) - U(:,2)))/norm(n12);

      R = [U(:,2) - U(:,4), U(:,1) - U(:,3), U(:,1) - U(:,5)];
      n13 = null(R');
      l13 = abs(n13'*(U(:,2) - U(:,1)))/norm(n13);

      R = [U(:,2) - U(:,5), U(:,1) - U(:,3), U(:,1) - U(:,4)];
      n14 = null(R');
      l14 = abs(n14'*(U(:,2) - U(:,1)))/norm(n14);

      R = [U(:,3) - U(:,5), U(:,1) - U(:,2), U(:,1) - U(:,4)];
      n15 = null(R');
      l15 = abs(n15'*(U(:,3) - U(:,1)))/norm(n15);

      [l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12 l13 l14 l15]

      Lmin = min([l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12 l13 l14 l15]);

    otherwise

      disp('Only 2, 3 and 4-dimensional input spaces are supported')

end

