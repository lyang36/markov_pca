%%get the next state. P: probability row
function s = getState(P)
        cumulative_distribution = cumsum(P);
        r = rand();
        s = find(cumulative_distribution > r,1);
end