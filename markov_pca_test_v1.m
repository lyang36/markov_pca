k1 = 50;
n1 = 2;
n = n1 * k1;
%deg = k1 - 1;

%generating orthostochastic matrix
%random regular graph
%p = full(createRandRegGraph(n, deg));
p = createRandClusterGraph(n1, k1);
deg = (diag(sum(p)));
p = deg\p;
%p = p / (deg);


transition_probabilities = p; 


starting_value = 1; 
chain_length = 10000;
chain = zeros(1, chain_length);
chain(1) =  starting_value;

skip_length = 1;

k = 2;

%jones = eye(n);
[u0, d0] = eig(p);
zt = u0(:, (k+1):n);
%diag(d0)
%u(:, 1)
%%%
eta0 = 0.01;
eta = eta0;
standard_base = eye(n);

s0 = standard_base(:, 1:k);
angle = zeros(1, chain_length);

iters = 0;
%last_time = 1;
this_time = 2;
for i = 2:skip_length:chain_length
    
    update_s = zeros(n, k);
    for j = 1:skip_length
        this_step_distribution = transition_probabilities(chain(i-1),:);
        cumulative_distribution = cumsum(this_step_distribution);
        r = rand();
        chain(this_time) = find(cumulative_distribution > r,1);
        update_s = update_s + eta * standard_base(:, chain(this_time - 1)) * s0(chain(this_time), :);
        this_time = this_time + 1;
    end
   
    %run sgd on the markov chain
    %s0(1)
    
    %if(mod(i, skip_length) == 0)
    iters = iters + 1;

    s0 = s0 + eta * standard_base(:, chain(i-1)) * s0(chain(i), :);
    %s0 = s0 + eta * p * s0;
    %size(s0)
    [us, rs] = qr(s0);
    s0(:, 1:k) = us(:, 1:k);
    ags = abs(s0' * u0(:, 1:k));
    %for j=1:k
    %    angle(j, i) = ags(j);
    %end
    angle(iters) = sum(sum((zt' * s0).^2));
    %angle(i) = abs(s0' * u0(:, 1:k));

    if (100 / iters < eta0)
        eta = 100/i;
    end
    %end
end
%chain
figure
hold on
d0(1)
plot(angle(1: iters))
% for i=1:k
%     plot(angle(i, :))
% end