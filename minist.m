tic;
images = loadMNISTImages('train-images-idx3-ubyte');
labels = loadMNISTLabels('train-labels-idx1-ubyte');
%Dig3 = images(:,labels==3);
%Dig8 = images(:,labels==8);
%[m,k] = size(Dig3);
%dig3 = Dig3(:,1:k);
%Dig3 = (Dig3'-mean(Dig3'))';
%dig8 = (dig8'-mean(dig8'))';
images=images(:,labels==3| labels==4| labels==5| labels==9);
[m,k] = size(images);
Cov = (1/(k-1))*images*images';
S = Cov(1:m/2, m/2+1:m);
[U,V] = eig([zeros(m/2,m/2),S;S',zeros(m/2,m/2)]);
% x0 = randn(m/2,1);
% y0 = randn(m/2,1);
% x0 = x0/norm(x0);
% y0 = y0/norm(y0);

%= Dig3(:,randperm(length(1:k)));
%     for l=1:5
%         images1=images(:,randperm(length(1:k)));
%            x1=x0;
%            y1=y0;
%         for i = 1:k
%          eta = 1e-4;
%           x = images1(1:m/2,i);
%           y = images1(m/2+1:m,i);
%           x2 = x1 + eta*x*y'*y1;
%           y2 = y1 + eta*y*x'*x1;
%           x1 = x2/norm(x2);
%            y1 = y2/norm(y2);
%   %  W(i,:) = 1/sqrt(2)*(U'*[x1;y1]);
%         W(l,i) = 1-abs(1/sqrt(2)*(U(:,m)'*[x1;y1]));
%         end;
%     end;
%     
%     Q=mean(W);
% 
%    plot(Q);
%    tsgd_1_4=toc;
for l=1:5
images1=images(:,randperm(length(1:k)));
  x1=x0;
  y1=y0;
  for i = 1:k
    eta = 0.05/i;
    x = images1(1:m/2,i);
    y = images1(m/2+1:m,i);
%     z2 = z1 + eta*[zeros(784,784),x*y';y*x',zeros(784,784)]*z1;
%     z1 = z2/norm(z2);
    x2 = x1 + eta*x*y'*y1;
    y2 = y1 + eta*y*x'*x1;
    x1 = x2/norm(x2);
    y1 = y2/norm(y2);
%    W(i,:) = 1/sqrt(2)*(U'*[x1;y1]);
     W(l,i) = 1-abs(1/sqrt(2)*(U(:,m)'*[x1;y1]));
  end
end;
   Q=mean(W);
 plot(Q);
 t_i=toc;
 
 
 