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
 [U0,V0] = eig(Z);
%% Random start or fixed start.
 %z0 = randn(6,1);
  z0=U0(:,5);
%% high dimension
 %z0 = randn(100,1);
%   z0 = U0(:,51)
z0 = z0/norm(z0);

R = chol(SS);
%% simulation loop
for j = 1:100
%% initializion
x1 = z0(1:3);
x1 = x1/norm(x1);
y1 = z0(4:6);
y1 = y1/norm(y1);
z1 = [x1;y1]/sqrt(2);
%% high dimension
% x1 = z0(1:50);
% x1 = x1/norm(x1);
% y1 = z0(51:100);
% y1 = y1/norm(y1);
% z1 = [x1;y1];
%% sgd loop

%% using the fixed stepsize
 for i = 1:n
    eta = 5e-5;
    z = R'*randn(6,1);
    x = z(1:3);
    y = z(4:6);
%     z2 = z1 + eta*[zeros(3,3),x*y';y*x',zeros(3,3)]*z1;
%     z1 = z2/norm(z2);
% %% for projection
%     x2 = x1 + eta*x*y'*y1;
%     y2 = y1 + eta*y*x'*x1;
%     x1 = x2/norm(x2);
%     y1 = y2/norm(y2);
%     W(i,:) = 1/sqrt(2)*(U0'*[x1;y1]);
%% for forest
z2 = z1+eta*([zeros(3,3),x*y';y*x',zeros(3,3)]*z1-z1'*[zeros(3,3),x*y';y*x',zeros(3,3)]*z1*z1);
z1 = z2;
    W(i,:) = (U0'*z1);
 end
 % plot each simulation;
    figure(1)
    plot(W(:,6));hold on;
    figure(2)
    plot(W(:,5));hold on;
    figure(3)
    plot(W(:,4));hold on;
    figure(4)
    plot(W(:,3));hold on;
    figure(5)
    plot(W(:,2));hold on;
    figure(6)
    plot(W(:,1));hold on;
    i=0;
   for k=[10,100,1000,50000,100000,150000,200000]
       KDE(j,6*i+1:6*i+6)=W(k,:)';
       i=i+1;
   end
   j
end;

%  for i=1:3
%  figure(4)
%  ksdensity(KDE(:,6*i-2));hold on;
%  figure(5)
%  ksdensity(KDE(:,6*i-1));hold on;
%  figure(6)
%  ksdensity(KDE(:,6*i));hold on;
%  end
%  for i=4:8
%  figure(7)
%  ksdensity(KDE(:,6*i-2));hold on;
%  figure(8)
%  ksdensity(KDE(:,6*i-1));hold on;
%  figure(9)
%  ksdensity(KDE(:,6*i));hold on;
%  end

