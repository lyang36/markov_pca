clear;
n = 2e5 ;
%% artificial or random generate matrix
%Sab = diag([0.3,0.26]); %for phase 1
%  Sab=diag(2*(rand(3,1)-0.5));
%  Sa = diag(diag(ones(3)));
%  Sb = diag(diag(ones(3)));
 
Sab = diag([4,2,0.5]);
Sa = [6,2,1;2,6,2;1,2,6];
Sb = Sa;
S = [Sa,Sab;Sab',Sb];
 
%% high dimension
% Sab = diag(99:-2:1)
% for i=1:50
%     for j=1:50
%         Sa(i,j)=300*2^(-abs(i-j));
%     end;
% end;
% Sb = Sa;
% S = [Sa,Sab;Sab',Sb];
%% Enlarge the singular gap.
 A = orth(randn(3));
 B = orth(randn(3));
%% high dimension.
%   A = orth(randn(50));
%   B = orth(randn(50));

 SSab = A'*Sab*B;
 SSa = A'*Sa*A;
 SSb = B'*Sb*B;
 SS = [SSa,SSab;SSab',SSb];
% SS = S;
%% Create the final matrix
 SSxy = SS(1:3,4:6);
 Z = [zeros(3,3),SSxy;SSxy',zeros(3,3)];
 
%% high dimension
%  SSxy = SS(1:50,51:100);
%  Z = [zeros(50,50),SSxy;SSxy',zeros(50,50)];
 [U,V] = eig(Z);
%% Random start or fixed start.
 %z0 = randn(6,1);
%  z0=U(:,5);
%% high dimension
 %z0 = randn(100,1);
%   z0 = U0(:,51)
%z0 = z0/norm(z0);

R = chol(SS);
%% simulation loop
%for j = 1:100
%% initializion
x0 = randn(3,1);
x0 = x0/norm(x0);
y0 = randn(3,1);
y0 = y0/norm(y0);
M0 = x0*y0';
%x0*y0';
M1 = M0;
% x1 = x0;
% y1 = y0;
  for i = 1:n
          eta = 0.05/i;
          z = R'*randn(6,1);
          x = z(1:3);
          y = z(4:6);
          M2 = M1+eta*x*y';
          [U0,S,V0] = svd(M2);
          v = diag(S);
          x = proj_L1_Linf(v,1);
          V1 = diag(x);
          M1 = U0*V1*V0';
  %  W(i,:) = 1/sqrt(2)*(U'*[x1;y1]);
          W(i) = 1-abs(U0(:,1)'*1/sqrt(2)*U(1:3,6)+V0(:,1)'*1/sqrt(2)*U(4:6,6));
  end;
  figure(3);
  plot(W);
toc;
