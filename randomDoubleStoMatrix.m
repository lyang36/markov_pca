%% generating a random stochastic matrix
%% each row is summed to 1
%% method: first randomly generate a matrix
%% then repeatedly make it row-stochastic and column stocastic

%% only stochastic matrix is produced if d1 ~= d2

function m = randomDoubleStoMatrix(d1, d2)
    m = randn([d1, d2]).^2;
    %eps = 1e-9;
    num_iter = 50;
    %while (sum(abs(sum(m, 1) - 1)> eps) >= 1) || (sum(abs(sum(m, 2) - 1) > eps) >= 1)
    for i=1:num_iter
        m = m/diag(sum(m, 1));
        %diag(sum(m, 1))
        m = diag(sum(m,2))\m;
        %diag(sum(m, 2))
    end
end