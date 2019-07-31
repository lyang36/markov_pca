clear;
n = 1e5;

Sab = diag([4,2,0.5]);
Sa = [6,2,1;2,6,2;1,2,6];
Sb = Sa;
S = [Sa,Sab;Sab',Sb];

A = orth(randn(3));
B = orth(randn(3));

 SSab = A'*Sab*B;
 SSa = A'*Sa*A;
 SSb = B'*Sb*B;
 SS = [SSa,SSab;SSab',SSb];

SSxy = SS(1:3,4:6);
Z = [zeros(3,3),SSxy;SSxy',zeros(3,3)];

 [U0,V0] = eig(Z);
 
   z0=U0(:,1)+0.0001*ones(6,1);
   z1 = z0/norm(z0);
   for i=1:100000
         eta = 0.00001;
      z1 = z1+ 2*eta*Z*z1;
      z1 = z1/norm(z1);
      W(i,:) = U0'*z1;
   end;
   figure(1)
   plot(W(:,1));
   figure(2)
   plot(W(:,2));
   figure(3)
   plot(W(:,3));
   figure(4)
   plot(W(:,4));
   figure(5)
   plot(W(:,5));
   figure(6)
   plot(W(:,6));
   