 d = 20; %dimension
 rk = 4; %rank
 chain_length = 50000; % number of samples
 isDrawSVD = 1;

%generating probability matrix
%random regular graph
p1 = randomStochasticMatrix(d, rk);
p2 = randomStochasticMatrix(rk, rk);
p3 = randomStochasticMatrix(rk, d);
p = p1 * p2 * p3;


transition_probabilities = p; 

% get the stationary distribution
[u, sigma] = eigs(p', 1);
u = u / sum(u);
stM = diag(1./u);
stM1 = diag(u);
p1 = stM1 * p;


% svd of p
[u0, sigma, v0] = svd(p);
[u1, sigma, v1] = svd(p1);

starting_value = 1; 
chain = zeros(1, chain_length);
chain(1) =  starting_value;

m0 = zeros(d);

values = zeros(1, chain_length);
values(1) = sum(sum((m0 - p).^2));
values1 = zeros(1, chain_length);
values1(1) = values(1);
%chain_matrices = zeros(d, d, chain_length);

for i=2:chain_length
    this_step_distribution = transition_probabilities(chain(i-1),:);
    cumulative_distribution = cumsum(this_step_distribution);
    r = rand();
    chain(i) = find(cumulative_distribution > r,1);
    
    m0(chain(i-1), chain(i)) = m0(chain(i-1), chain(i)) + 1;
    m1 = stM * m0 / i;
    
    if(isDrawSVD)
        %chain_matrices(:,:, i) = m1;
        [us, sigmas, vs] = svd(m1);
        angs = us(:, 1:rk)' * u0(:, 1:rk);

        [us1, sigmas1, vs1] = svd(m0);
        angs1 = us1(:, 1:rk)' * u1(:, 1:rk);
        values(i) = rk - trace(angs' * angs);
        values1(i) = rk -trace(angs1' * angs1);

        %values(i) = sqrt(sum(sum((m1 - p).^2)));
    end
end
%chain
if(isDrawSVD)
    figure
    plot(1:chain_length, values)
    hold on
    plot(1:chain_length, values1)
    % 
    legend('(sin Theta)^2: svd(P)', '(sin Theta)^2: svd(M*P)')
    title(sprintf('rank = %d, dimension=%d', rk, d))
    hold off
end