%%% run msg for num_repeats times, and output the average
%%% using fix step size

deg = 3;
n = n1;
transition_probabilities = p0; 
starting_value = 1; 
skip_length = 2;


chain_length = 30000;
xi=1;

eta0 = 0.0075;
%eta = 0.04;
eta = 0.0075

k = deg;

num_repeats = 20;
angle_msg = zeros(num_repeats, chain_length);
time_msg = zeros(num_repeats, chain_length);
iters_msg = zeros(2*n, k, num_repeats, chain_length);

for rts = 1:num_repeats
    x1 = zeros(n);
    t1 = now;
    for i=1:chain_length
        for rskip = 1:skip_length
            xi0 = xi;
            this_step_distribution = transition_probabilities(xi0,:);
            cumulative_distribution = cumsum(this_step_distribution);
            r = rand();
            xi = find(cumulative_distribution > r,1);
        end

        online_sample = zeros(n);
        online_sample(xi0, xi) = 1;

        %eta = 0.0075*log(i) / i;
    %     if(mod(i, 300) == 1 && i > 300 && isDiminishing_step == 1)
    %         eta =  300 * eta0 / i / 2 ;
    %         disp(eta)
    %     end
        if(mod(i, 6000) == 1)
            eta = 6000 * eta0 / (6000 + i);
        end
        %%%%
        % use online_sample to do the MSG
        x2 = x1 + eta * online_sample;
        [U0, S, V0] = svd(x2);
        v = diag(S);
        x = proj_L1_Linf(v, k);
        V1 = diag(x);
        x1 = U0 * V1 * V0';
        
        iters_msg(:,:, rts, i) = [U0(:, 1:k); V0(:, 1:k)];
        time_msg(rts, i) =  (now - t1) * 3600 * 24 * 1000; %milliseconds 
        %angle_msg(rts, i) = k - (norm(U0(:, 1:k)' * us1(:, 1:k), 'fro')^2 + norm(V0(:, 1:k)' * vs1(:, 1:k), 'fro')^2)/2;
    end
end
%semilogy(angle_msg)



%compute angles
for rpt = 1:num_repeats 
    for i=1:chain_length
        s0 = iters_msg(:,:, rpt, i);
        U0 = s0(1:n, :);
        V0 = s0((n+1):(2*n), :);
        angle_msg(rpt, i) = k - (norm(U0(:, 1:k)' * us1(:, 1:k), 'fro')^2 + norm(V0(:, 1:k)' * vs1(:, 1:k), 'fro')^2)/2;
    end
end

figure
hold on
for i=1:num_repeats
    semilogy(time_msg(i, :), angle_msg(i, :));
end

semilogy(sum(time_msg, 1)/num_repeats, sum(angle_msg, 1)/num_repeats)
