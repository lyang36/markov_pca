chain_length = 5000; % number of samples

isDrawSVD = false;
%d = 40; %dimension
%rk = 1; %rank


% fix sample size, fix dimension 
d = 100;

erros_rank_error = zeros(1, d);
for rk = 1:d
   rk
   run_markov_sample
   erros_rank_error(rk) = sqrt(sum(sum((m1 - p).^2)));
end

figure
plot(1:d, erros_rank_error);
%legend('Fix sample size 6000, dimension = 100, grow rank')
legend(sprintf('Fix sample size %d, fix dimension = %d, grow rank', chain_length, d))
xlabel('rank')
ylabel('||M-P||_F')


% fix sample size, fix rank
rk = 10;

erros_dim_error = zeros(1, d-rk);
for d = rk : 100
    d
   run_markov_sample
   erros_dim_error(d - rk + 1) = sqrt(sum(sum((m1 - p).^2)));
end

figure
plot(rk : 100, erros_dim_error);
%legend('Fix sample size 6000, fix rank = 10, grow dimension')
legend(sprintf('Fix sample size %d, fix rank = %d, grow dimension', chain_length, rk))
xlabel('dimension')
ylabel('||M-P||_F')

% fix sample size, fix dimension
% do a svd for the matrix, and then measure the angle
d = 100;

sin_theta_sqr_error = zeros(1, d);
for rk = 1:d
   rk
   run_markov_sample
   [us, sigmas, vs] = svd(m0);
   angs = us(:, 1:rk)' * u1(:, 1:rk);
   sin_theta_sqr_error(rk) = rk - trace(angs' * angs);
end

figure
plot(1:d, sin_theta_sqr_error);
legend(sprintf('Fix sample size %d, dimension = %d, grow rank', chain_length, d))
xlabel('rank')
ylabel('(sin Theta(r))^2')

% 
