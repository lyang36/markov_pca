%generate graph for compare convex, non-convex
% comment the eta setting before runing this script

markov_sgd_diminishing_eta
markov_pca_msg_diminishing_eta

% 
% eta = 0.05
% markov_pca_test_input_markov_sample
% angles_non_convex_005 = angle_all(1, :);
% markov_pca_msg
% angles_msg_005 = angle_msg;
% 
% 
% 
% eta = 0.
% markov_pca_test_input_markov_sample
% angles_non_convex_001 = angle_all(1, :);
% markov_pca_msg
% angles_msg_001 = angle_msg;




figure
z = isfinite(angle_all(:, length(angle_all(1,:))))
angles_non_convex_dim_01 = (sum(angle_all(z, :), 1)/(repeats-1));
angles_msg_dim_01 = sum(angle_msg, 1)/num_repeats;

semilogy(sum(time_all(z, :), 1)/(repeats-1), angles_non_convex_dim_01, '-')
hold on
% semilogy(angles_non_convex_005, '-.')
% semilogy(angles_non_convex_001, '--')
semilogy(sum(time_msg, 1)/num_repeats, angles_msg_dim_01, '-')
% semilogy(angles_msg_005, '-.')
% semilogy(angles_msg_001, '--')
%legend('SGD \eta=0.1', 'SGD \eta=0.05', 'SGD \eta=0.01','MSG \eta=0.1', 'MSG \eta=0.05', 'MSG \eta=0.01')
legend('SGD', 'MSG')
xlabel('Milliseconds')
ylabel('$\|\sin\Theta\|_F^2$', 'interpreter', 'latex')
%axis([0, 1e4, 0.05, 2])

save run_comparison_dim_step p0 pt n1 deg angles_non_convex_01 angles_msg_01 time_all time_msg repeats num_repeats