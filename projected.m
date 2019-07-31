function [S,flag]=projected(eigvalue,multi,d,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Introduction
% programmed by Zhehui Chen(Gatech)
% email: zhchen@gatech.edu
% last modified Feb. 20th 2017
% implenment an effiencient algorithm to do matrix projection
%% Input
% eigvalue : the distinct eigen-value.
% multi : multiplicities of the eigen-value.
% d : dimension of eigen-vector.
% k : rank of matrix.
%% Output
% sigma : the projection eigen-value.
% flag : [ 0 means error ]
%        [ 1 means exist sigma]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,e] = size(eigvalue);
[sigma,I] = sort(eigvalue);
kappa = multi(I);
flag = 1;
i = 1;
j = 1;
s1 = 0;
s2 = 0;
c1 = 0;
c2 = 0;
while( i <= n)
    if(i<j)
        S = (k-(s2-s1)-(d-c2))/(c2-c1);
        b = ((sigma(i)+S>=0) && (sigma(j-1)+S<=1) && ((i<=1) || (sigma(i-1)+S<=0)) && ((j>=n)) || (sigma(j+1)>=1));
        if b;
            return ;
        end;
        if ((j<=n) && (sigma(j)-sigma(i) <=1))
           s2 = s2 +kappa(j)*sigma(j);
           c2 = c2 +kappa(j);j=j+1;
        else
           s1 = s1 +kappa(i)*sigma(i);
           c1 = c1 +kappa(i);i=i+1;
        end;   
    end;
end;
flag = 0;
return;