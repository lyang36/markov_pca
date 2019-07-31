n = 15;
d = 10;

A = randn(n,d);
S = A'*A/100;

[U0,V0] = eig(S);
v = diag(V0);
x = proj_L1_Linf(v,1);
V1 = diag(x);

Q = U0*V1*U0';