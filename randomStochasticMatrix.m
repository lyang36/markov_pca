%% generating a random stochastic matrix
%% each row is summed to 1

function m = randomStochasticMatrix(d1, d2)
    m = zeros(d1, d2);
    for i = 1 : d1
        u = randn(1, d2);
        u = u / norm(u);
        u = u.^2;
        m(i, :) = u;
    end
end