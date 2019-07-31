n = 80;
deg = 10;

%generating orthostochastic matrix
%random regular graph
p = full(createRandRegGraph(n, deg));
p = p / (deg);


transition_probabilities = p; 


starting_value = 1; 
chain_length = 50000;
chain = zeros(1, chain_length);
chain(1) =  starting_value;

skip_length = 1;

k = 2;

%jones = eye(n);
[u0, d0] = eigs(p, n, 'sr');

zt = u0(:, (k+1):n);
%diag(d0)
%u(:, 1)
%%%
eta0 = 0.1;
eta = eta0;
standard_base = eye(n);

s0 = standard_base(:, 1:k);
angle = zeros(1, chain_length);
iters = 0;
for i=2:chain_length
    this_step_distribution = transition_probabilities(chain(i-1),:);
    cumulative_distribution = cumsum(this_step_distribution);
    r = rand();
    chain(i) = find(cumulative_distribution > r,1);
    

    
    %run sgd on the markov chain
    %s0(1)
    
    if(mod(i, skip_length) == 0)
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
        
        if (2000 / iters < eta0)
            eta = 2000/i;
        end
    end
end
%chain
figure
hold on
d0(1)
plot(angle(1: iters))
% for i=1:k
%     plot(angle(i, :))
% end