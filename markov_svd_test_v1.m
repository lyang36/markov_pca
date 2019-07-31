%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%this code solves the problem of 
%%svd the transition matrix

d = 60; %dimension
rk = 4; %rank

%generating probability matrix
%random regular graph
p1 = randomStochasticMatrix(d, rk);
p2 = randomStochasticMatrix(rk, rk);
p3 = randomStochasticMatrix(rk, d);
p = p1 * p2 * p3;

n = d;

%random regular graph
%p = full(createRandRegGraph(n, deg));
% p = createRandClusterGraph(n1, k1);
% deg = (diag(sum(p)));
% p = deg\p;
%p = p / (deg);



transition_probabilities = p; 


starting_value = 1; 
chain_length = d * 6000;
chain = zeros(1, chain_length);
chain(1) =  starting_value;

skip_length = 1;

k = 1;

%jones = eye(n);
[v11, d11] = eigs(p',1);
stp = diag(abs(v11) / sum(abs(v11)));
[u0, d0, v0] = svd(stp * p);

solutions = [u0(:, 1:k); v0(:, 1:k)];
zt = [u0(:, (k+1):n); v0(:, (k+1):n)] / sqrt(2);
dpmat = [zeros(n), p; p', zeros(n)];

%%%
eta0 = 0.001;
eta = eta0;
standard_base = eye(n);

s0 = [standard_base(:, 1:k); standard_base(:, 1:k)]/sqrt(2);
angle = zeros(1, chain_length);

iters = 0;
%last_time = 1;
this_time = 2;
current_state = starting_value;

for i = 2:skip_length:chain_length
   
    %run sgd on the markov chain
    %independent
    state1 = getState(diag(stp)');
    %nonindependent
    %state1 = getMarkovState(p, current_state);
    current_state = getMarkovState(p, state1);
    iters = iters + 1;

    %xmat = standard_base(:, state1) * standard_base(:, current_state)';
    %zmat = [zeros(n), xmat; xmat', zeros(n)];
    
    %ds0 = s0 + eta * zmat * s0;
    ds0 = s0 + eta * [standard_base(:, state1) * s0(n + current_state, :); ...
                            standard_base(:, current_state) * s0(state1, :)];
    s0 = s0 + ds0 / norm(ds0);

    [us, rs] = qr(s0);
    s0(:, 1:k) = us(:, 1:k);
    ags = abs(s0' * zt(:, 1:k));

    angle(iters) = sum(sum((zt' * s0).^2));

    %if (100 / iters < eta0)
    %    eta = 100/i;
    %end
    %end
end
%chain
figure
hold on
d0(1)
plot(angle(1: iters))