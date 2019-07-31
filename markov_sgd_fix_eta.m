%%run markov SGD for number_repeats times, using fix step size
%% output the average convergence curve,
%% used to compare the convergence with MSG
deg = 3;
n = n1 * 2;
ll = abs(ll)
p = [zeros(n/2), pt; pt', zeros(n/2)];
transition_probabilities = p0; 
starting_value = 1; 
skip_length = 2;

k = deg;

zt = ([us1(:, (k+1):(n/2)); vs1(:, (k+1):(n/2))])/sqrt(2);
z0 = orth([us1(:, 1:k); vs1(:, 1:k)]);

eta = 0.04;
errthres = eta;
chain_length = 200000;
standard_base = eye(n / 2);

%s0 = [standard_base(:, 1:k); standard_base(:, 1:k)] / sqrt(2);
repeats = 20;
angle_all = zeros(repeats, chain_length);
time_all = zeros(repeats, chain_length);
iters_all = zeros(n, k, repeats, chain_length);

initial_all = zeros(repeats, 1);


for rpt = 1:repeats 
    s0 = orth(randn(n, k));
    %angle = zeros(1, chain_length);
    iters = 0;
    xi = 1;

    initial_all(rpt, 1) = norm(zt' * s0 / (s0' * z0), 'fro')^2;

    emp_u = zeros(n/2, 1);
    
    t1 = now;
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
        
        iters_all(:,:, rpt, i) = s0;
        
        time_all(rpt, i) = (now - t1) * 3600 * 24 * 1000; %milliseconds 
        %[sz1, ~] = qr(s0, 0);
        %errang =  norm(zt' * sz1, 'fro')^2;
        
        %if(iters > 10)
        %    eta = 1000*0.1/(990 + iters)
        %end
        %if(iters > 300)
        %    eta = 0.12;
        %end
        
        
        %if(norm(s0) > k^2)
        %    [s0, ~] = qr(s0, 0);
        %end
        
        %angle(iters) = errang;
    end
    %angle_all(rpt, :) = angle;
end


%compute angles
for rpt = 1:repeats 
    angle = zeros(1, chain_length);
    for i=1:chain_length
        s0 = iters_all(:,:, rpt, i);
        [sz1, ~] = qr(s0, 0);
        errang =  norm(zt' * sz1, 'fro')^2;
        angle(i) = errang;
    end
    angle_all(rpt, :) = angle;
end
%%%
%time complexity
semilogy(sum(time_all, 1)/repeats, sum(angle_all, 1)/repeats)

