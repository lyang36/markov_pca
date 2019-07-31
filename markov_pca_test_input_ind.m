%%test the algorithm in the notes


n = 10;
deg = 3;

%generating orthostochastic matrix
%random regular graph
%p = full(createRandRegGraph(n, deg));
%p = p / (deg);
%random regular graph
p1 = randomStochasticMatrix(n/2, deg);
p2 = randomStochasticMatrix(deg, deg);
p3 = randomStochasticMatrix(deg, n/2);
p0 = p1 * p2 * p3;
p = [zeros(n/2), p0; p0', zeros(n/2)];


transition_probabilities = p0; 

ll = eigs(p, deg+1, 'la')

starting_value = 1; 
%chain = zeros(1, chain_length);
%chain(1) =  starting_value;

skip_length = 1;

k = 3;

%jones = eye(n);
[u0, d0] = eigs(p, n, 'la');

zt = u0(:, (k+1):n);
%diag(d0)
%u(:, 1)

% 500000;%

%%%
eta0 = 0.1;
eta = ( ll(k) - ll(k+1) ) / 30;

chain_length = round(10*log(n * k)/ eta / (ll(k) - ll(k+1))^2)


standard_base = eye(n / 2);

%s0 = [standard_base(:, 1:k); standard_base(:, 1:k)] / sqrt(2);
s0 = orth(randn(n, k));
angle = zeros(1, chain_length);
angle_emp = zeros(1, chain_length);

iters = 0;

empirical_matrix = zeros(n, n);

xi = 1;
for i=1:chain_length
    xi0 = randi(n/2);   %independent samples
    %xi0 = xi;
    this_step_distribution = transition_probabilities(xi0,:);
    cumulative_distribution = cumsum(this_step_distribution);
    r = rand();
    
    %chain(i) = find(cumulative_distribution > r,1);
    xi = find(cumulative_distribution > r,1);
   
    %run sgd on the markov chain
    %s0(1)
    
    iters = iters + 1;
    %s0 = s0 + eta * standard_base(:, xi0) * s0(xi, :);
    exp_p_0 = standard_base(:, xi0) * standard_base(:, xi)';
    exp_p = [zeros(n/2), exp_p_0; exp_p_0', zeros(n/2)];
    empirical_matrix = empirical_matrix + exp_p;
    
    s0 = s0 + eta * ( exp_p * s0 - s0 * (s0' * exp_p * s0));
    
    %eta = eta0/(10+ iters);
    
    %s0 = s0 + eta * p * s0;
    %size(s0)
    
    %[us, rs] = qr(s0);
    %s0(:, 1:k) = us(:, 1:k);
    %ags = abs(s0' * u0(:, 1:k));

    %[s1, ~] = eigs(empirical_matrix, k, 'la');
    angle(iters) = norm((zt' * s0).^2, 'fro')^2;
    
    if(iters > 1)
        angle_emp(iters) = angle_emp(iters - 1);
    end
    if(mod(i, 1000) == 0)
        [us, ds] = eigs(empirical_matrix / iters, k, 'la');
        %s0(:, 1:k) = us(:, 1:k);
        %ags = abs(s0' * u0(:, 1:k));
        angle_emp(iters) = norm((zt' * us(:, 1:k)).^2, 'fro')^2;
    end
    %angle_emp(iters) = sum(sum((zt' * s1).^2));
    %angle(i) = abs(s0' * u0(:, 1:k));

    %if (500 / iters < 0.01)
    %    eta = 500/i;
    %end
    
    if(mod(iters, 10000) == 0)
        plot(angle(1: iters));
        drawnow
    end
    
end
%chain
figure
hold on
d0(1)
lx =(1: iters);
plot(log10(angle(1: iters)))
plot(log10(angle_emp(1: iters)))

%plot((angle_emp(1: iters)))
plot(log10(angle(1)*exp(-(ll(k)-ll(k+1)) / (1.5 * n) * eta * lx)))
%str = '$$N=\tilde{O}\left(\frac{d}{\epsilon(\lambda_{k}-\lambda_{k+1})}\right), \epsilon=1/5$$';
%text(0.25,2.5,str,'Interpreter','latex')
%title(str,'Interpreter','latex')
%text(chain_length/3,1,'Theoretic bound: $N=\tilde{O}\left(\frac{d}{\epsilon(\lambda_{k}-\lambda_{k+1})^{3/2}}\right)$','Interpreter','latex')

legend('Algorithm', 'empirical estimation', 'theory')
norm(us(:, 1:k)'*u0(:, 1:k),'fro')^2
norm(orth(s0)'*u0(:, 1:k),'fro')^2
ll(k) - ll(k+1)
