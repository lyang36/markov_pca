%%% this code is writen for Zhehui to run MSG
%%% inside the loop, online_sample is a n*n matrix

deg = 3;
n = 4 * deg ;
transition_probabilities = p0; 
starting_value = 1; 
skip_length = 2;


chain_length = 7000;
xi=1;

%eta0 = 0.1;
eta = 0.1;

x1 = zeros(n);
angle_msg = [];
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
    
    %eta = 0.1*log(i) / i;
%     if(mod(i, 300) == 1 && i > 300 && isDiminishing_step == 1)
%         eta =  300 * eta0 / i / 2 ;
%         disp(eta)
%     end
    %%%%
    % use online_sample to do the MSG
    x2 = x1 + eta * online_sample;
    [U0, S, V0] = svd(x2);
    v = diag(S);
    x = proj_L1_Linf(v, 3);
    V1 = diag(x);
    x1 = U0 * V1 * V0';
    angle_msg = [angle_msg, deg - norm(U0(:, 1:3)' * us1(:, 1:3), 'fro')^2];
end
semilogy(angle_msg)
