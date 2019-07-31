%this document illustrate the Euler's method of solving 
%the ODE

n = 80;
deg = 20;

%generating orthostochastic matrix
%random regular graph
p = full(createRandRegGraph(n, deg));
p = p / (deg);


transition_probabilities = p; 


chain_length = 5000;

%this is the subspace dimension
k = 2;

%jones = eye(n);
[u0, d0, v0] = svd(p);

solutions = [u0(:, 1:k); v0(:, 1:k)];
zt = [u0(:, (k+1):n); v0(:, (k+1):n)];
dpmat = [zeros(n), p; p', zeros(n)];

%starting step size
eta0 = 0.01;
%step size maybe further decreasing with the time steps
eta = eta0;
standard_base = eye(n);

s0 = [standard_base(:, 1:k); standard_base(:, 1:k)];
angle = zeros(1, chain_length);

iters = 0;
for i=2:chain_length
    iters = iters + 1;

    s0 = s0 + eta * ( dpmat * s0 - 1/2 * trace(s0' * dpmat * s0)*s0);
    
    angle(iters) = sum(sum((zt' * s0).^2));

    if (2000 / iters < eta0)
        eta = 2000/i;
    end
end
%chain
figure
hold on
plot(angle(1: iters))
