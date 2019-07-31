tic;
%clear;
% images = loadMNISTImages('train-images-idx3-ubyte');
% labels = loadMNISTLabels('train-labels-idx1-ubyte');
% images1=images(:,labels==3| labels==4| labels==5| labels==9);
% [m,k] = size(images1);
% 
% Cov = (1/(k-1))*images1*images1';
% S = Cov(1:m/2, m/2+1:m);
% [U,V] = eig([zeros(m/2,m/2),S;S',zeros(m/2,m/2)]);

%x0 = randn(m/2,1);
% x0 = [1;zeros(m/2-1,1)];
% %y0 = randn(m/2,1);
% y0 = [1;zeros(m/2-1,1)];
% %y0/norm(y0);
M0 = x0*y0';
%x0*y0';
for j=1:5
    
images1=images(:,randperm(length(1:k)));
  M1 = M0;
% x1 = x0;
% y1 = y0;
  for i = 1:k
          eta = 0.05*log(k)/k;
          x = images1(1:m/2,i);
          y = images1(m/2+1:m,i);
          M2 = M1+eta*x*y';
          [U0,S,V0] = svd(M2);
          v = diag(S);
          x = proj_L1_Linf(v,1);
           V1 = diag(x);
           M1 = U0*V1*V0';
  %  W(i,:) = 1/sqrt(2)*(U'*[x1;y1]);
           W(j,i) = 1-abs(U0(:,1)'*1/sqrt(2)*U(1:m/2,m)+1/sqrt(2)*V0(:,1)'*U(m/2+1:m,m));
  end;
  j
end;
Q=mean(W);
  plot(Q);
                      tmsg2=toc;
