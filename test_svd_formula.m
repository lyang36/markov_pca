%%test the algorithm in the notes
%%test whether the SVD formula works
%%the result is that it may not work

n = 10;
deg = 4;

%generating orthostochastic matrix
%random regular graph
%p0 = full(createRandRegGraph(n/2, deg));
%p = p / (deg);
%random regular graph
p1 = randomStochasticMatrix(n/2, deg);
p2 = randomStochasticMatrix(deg, deg);
p3 = randomStochasticMatrix(deg, n/2);
p0 = p1 * p2 * p3;


transition_probabilities = p0; 

%ll = eigs(p, deg+1, 'la')

starting_value = 1; 
%chain = zeros(1, chain_length);
%chain(1) =  starting_value;

skip_length = 1;

k = 2;

%jones = eye(n);
[u0, d0, v0] = svd(p0);

zt = u0(:, (k+1):(n/2));
%diag(d0)
%u(:, 1)

chain_length = 60000%round(10*log(n * k)/ eta / (ll(k) - ll(k+1))^2)

%%%
eta0 = 1;
eta = 0.01;

standard_base = eye(n / 2);

su0 = orth(randn(n/2, k));%standard_base(:, 1:k);
sv0 = orth(randn(n/2, k))%standard_base(:, 1:k);

angle = zeros(1, chain_length);
angle_emp = zeros(1, chain_length);

iters = 0;

empirical_matrix = zeros(n, n);

for i=2:chain_length
    xi0 = randi(n/2);
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
    %exp_p = [zeros(n/2), exp_p_0; exp_p_0', zeros(n/2)];
    %empirical_matrix = empirical_matrix + exp_p;
    
    su0 = su0 + eta * ( exp_p_0 * sv0 - su0 * (su0' * exp_p_0 * sv0));
    sv0 = sv0 + eta * ( exp_p_0' * su0 - sv0 * (sv0' * exp_p_0' * su0));

    %s0 = s0 + eta * p * s0;
    %size(s0)
    
    %su0 = su0 + eta * ( exp_p_0 * sv0);
    %sv0 = sv0 + eta * ( exp_p_0' * su0);
    %su0 = su0 + eta * ( p0 * sv0);
    %sv0 = sv0 + eta * ( p0' * su0);
    [us0, ~] = qr(su0);
    [vs0, ~] = qr(sv0);
    %su0 = us0(:, 1:k);
    %sv0 = vs0(:, 1:k);
    us0 = us0(:, 1:k);
    vs0 = vs0(:, 1:k);
    su0 = su0 + eta * ( exp_p_0 * sv0 - us0 * (us0' * exp_p_0 * sv0));
    sv0 = sv0 + eta * ( exp_p_0' * su0 - vs0 * (vs0' * exp_p_0' * su0));

    %s0(:, 1:k) = us(:, 1:k);
    %ags = abs(su0' * u0(:, 1:k));

    %[s1, ~] = eigs(empirical_matrix, k, 'la');
    angle(iters) = trace(su0' * p0 * sv0);%norm([su0; sv0]' * [u0(:,(k+1):(n/2)); v0(:,(k+1):(n/2))], 'fro')^2;%trace(su0' * p0 * sv0); %
    %angle_emp(iters) = sum(sum((zt' * s1).^2));
    %angle(i) = abs(s0' * u0(:, 1:k));

    %if (500 / iters < eta0)
    %    eta = 500/i;
    %end
    
end
%chain
figure
hold on
d0(1)
lx =(1: iters);
hold off
plot((angle(1: iters)))
%plot((angle_emp(1: iters)))
%plot(angle(1)*exp(-(ll(k)-ll(k+1)) / (1.5 * n) * eta * lx))
%str = '$$N=\tilde{O}\left(\frac{d}{\epsilon(\lambda_{k}-\lambda_{k+1})}\right), \epsilon=1/5$$';
%text(0.25,2.5,str,'Interpreter','latex')
%title(str,'Interpreter','latex')
%text(chain_length/3,1,'Theoretic bound: $N=\tilde{O}\left(\frac{d}{\epsilon(\lambda_{k}-\lambda_{k+1})^{3/2}}\right)$','Interpreter','latex')
hold on
d00 = diag(d0)
plot(ones([1, iters]) * sum(d00(1:k)))
legend('Algorithm', 'empirical estimation', 'theory')
svd(p0)
svd(su0'*p0*sv0)
