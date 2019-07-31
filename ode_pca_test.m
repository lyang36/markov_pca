%this document illustrate the Euler's method of solving 
%the ODE

n = 80;
deg = 10;

%generating orthostochastic matrix
%random regular graph
%p = full(createRandRegGraph(n, deg));
%p = p / (deg);
p1 = randomStochasticMatrix(n/2, deg);
p2 = randomStochasticMatrix(deg, deg);
p3 = randomStochasticMatrix(deg, n/2);
p = p1 * p2 * p3;
p = [zeros(n/2), p; p', zeros(n/2)];


transition_probabilities = p; 


chain_length = 5000;

%this is the subspace dimension
k = 4;

%jones = eye(n);
[u0, d0] = eigs(p, n, 'sr');

zt = u0(:, (k+1):n);
zt0 = u0(:, 1:k);

%starting step size
eta0 = 0.01;
%step size maybe further decreasing with the time steps
eta = eta0;
standard_base = eye(n);

s0 = standard_base(:, 1:k);
s1 = standard_base(:, 1:k);
angle = zeros(1, chain_length);
angle1 = zeros(1, chain_length);


iters = 0;
for i=2:chain_length
    iters = iters + 1;

    s1 = s1 + eta * ( p * s1 - s1 * (s1' * p * s1));
    
    %oja
    s0 = s0 + eta * p * s0;
    [us, rs] = qr(s0);
    s0(:, 1:k) = us(:, 1:k);

    angle(iters) = sum(sum((zt0(:,1)' * s1).^2));
    angle1(iters) = sum(sum((zt' * s1).^2))/sum(sum((zt0' * s1).^2));

    %if (2000 / iters < eta0)
    %    eta = 2000/i;
    %end
end
%chain
hold off
plot(angle(1: iters))
hold on
plot(angle1(1: iters))

lambda0=d0(k-1,k-1)
lambda1=d0(k,k)
plot((1: iters), n / 4 *exp(-(1: iters)*eta*(lambda0-lambda1)))

