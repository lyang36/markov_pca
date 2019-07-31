%%test the algorithm in the notes
%%this code full depends on the random walk generated from the chain.


deg = 3;
n = 4 * deg * 2;

%%%%initialization%%%
%generating orthostochastic matrix
%random regular graph
%p = full(createRandRegGraph(n, deg));
%p = p / (deg);
%random regular graph
% p1 = randomStochasticMatrix(n/2, deg);
% p2 = randomStochasticMatrix(deg, deg);
% p3 = randomStochasticMatrix(deg, n/2);
% p0 = p1 * p2 * p3;
% pm = randperm(n/2);
% p0 = createRandLumpGraph(n/2, deg, 0.5);
% p0 = p0(pm, pm)
% 
% [u0, d0] = eig(p0');
% u0 = u0(:,1) / sum(u0(:,1));
% ll = eigs(diag(u0) * p0)
ll = abs(ll)

p = [zeros(n/2), pt; pt', zeros(n/2)];


transition_probabilities = p0; 


starting_value = 1; 
skip_length = 2;

k = deg;

zt = orth([us1(:, (k+1):(n/2)); vs1(:, (k+1):(n/2))]);
z0 = orth([us1(:, 1:k); vs1(:, 1:k)]);

%%%
%eta0 = 0.1;
%eta = 0.1;

errthres = eta;


chain_length = 7000;%round(sum(ll)*log(n * k)/ eta / (ll(k) - ll(k+1))^2 / 10)


standard_base = eye(n / 2);

%s0 = [standard_base(:, 1:k); standard_base(:, 1:k)] / sqrt(2);
repeats = 1;
angle_all = zeros(repeats, chain_length);
initial_all = zeros(repeats, 1);


for rpt = 1:repeats 
    s0 = orth(randn(n, k));
    angle = zeros(1, chain_length);
    angle_emp = zeros(1, chain_length);
    iters = 0;
    empirical_matrix = zeros(n, n);
    xi = 1;

    initial_all(rpt, 1) = norm(zt' * s0 / (s0' * z0), 'fro')^2;

    emp_u = zeros(n/2, 1);
    for i=1:chain_length
        %xi0 = randi(n/2);   %independent samples
        for rskip = 1:skip_length
            xi0 = xi;
            this_step_distribution = transition_probabilities(xi0,:);
            cumulative_distribution = cumsum(this_step_distribution);
            r = rand();
            %get the next sample
            xi = find(cumulative_distribution > r,1);
        end

        iters = iters + 1;
        emp_u(xi) = emp_u(xi) + 1;


        %%%sparse update
        pid = xi0;
        did = xi;
        n1 = n/2;
        w1 = (s0 * s0(pid, :)') * s0(n1 + did, :) + (s0 * s0(n1 + did, :)') * s0(pid, :);
        s0(pid, :) = s0(pid, :) + eta * s0(n1 + did, :);
        s0(n1 + did, :) = s0(n1 + did, :) + eta * s0(pid, :);
        s0 = s0 - eta * w1;
        empirical_matrix(pid, n1 + did) = empirical_matrix(pid, n1 + did) + 1;
        empirical_matrix(n1 + did, pid) = empirical_matrix(n1 + did, pid) + 1;


        
        %if(iters > 1)
        %    angle(iters) = angle(iters - 1);
        %    angle_emp(iters) = angle_emp(iters - 1);
        %end

        % compute the subspace angle
        %if(mod(i, 100) == 1)
        %    [us, ds] = eigs(empirical_matrix, k, 'la');
        %    angle_emp(iters) = norm((zt' * us(:, 1:k)), 'fro')^2;

            [sz1, ~] = qr(s0, 0);
            errang =  norm(zt' * sz1, 'fro')^2;
            angle(iters) = errang;
        %end

        %eta = 0.1*log(i) / i;
%         if(mod(i, 500) == 1 && i > 500 && isDiminishing_step == 1)
%             eta =  500 * eta0 / i / 2 ;
%             disp(eta)
%         end

%         if(mod(iters, 10000) == 0)
%             hold off
%             plot(angle(1: iters));
%             hold on
%             plot(angle_emp(1: iters));
%             grid on;
%             drawnow
%         end

    end
    angle_all(rpt, :) = angle;
end
%chain
figure
lx =(1: 10^5);
init_exp = median(initial_all);
for rpt = 1:repeats
    semilogy((angle_all(rpt, 1: chain_length)), '--')
    hold on
end
semilogy((init_exp*exp(-abs((ll(k)-ll(k+1)) / ll(1)) * eta * lx / 4)), '-k')
axis([1, 7e4, 0.002, 10])
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'--r');
h(2) = plot(NaN,NaN,'-k');
legend(h, 'red','blue','black');
legend(h, 'Algorithm', 'Thoeretic upper bound')%'empirical estimation', 'theory')
ylabel('$$\|\sin\Theta\|_F^2$$', 'Interpreter', 'latex')
xlabel('Iterations')

%group_learn = kmeans(diag(emp_u)\us_learn, 3)