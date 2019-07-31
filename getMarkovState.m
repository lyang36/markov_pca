%%get the next markov state. P: transition matrix
%%s0 -- current state

function s = getMarkovState(P, s0)
        this_step_distribution = P(s0,:);
        cumulative_distribution = cumsum(this_step_distribution);
        r = rand();
        s = find(cumulative_distribution > r,1);
end