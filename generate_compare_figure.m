load run_comparison_fix_step.mat

figure
subplot(2,1,1)
semilogy(sum(time_all, 1)/repeats, angles_non_convex_01, '--r')
hold on
semilogy(sum(time_msg, 1)/num_repeats, angles_msg_01, '-k')
%xlabel('Milliseconds')
ylabel('$\|\sin\Theta\|_F^2$', 'interpreter', 'latex')
axis([0, 5000, 0.1, 2])
set(gca,'position',[0.1 0.55 0.85 0.43],'units','normalized')
xticks([])
xlabel('')

load run_comparison_dim_step.mat
subplot(2,1,2)
z = isfinite(angle_all(:, length(angle_all(1,:))))
angles_non_convex_dim_01 = (sum(angle_all(z, :), 1)/(repeats-1));
angles_msg_dim_01 = sum(angle_msg, 1)/num_repeats;

semilogy(sum(time_all(z, :), 1)/(repeats-1), angles_non_convex_dim_01, '--r')
hold on
semilogy(sum(time_msg, 1)/num_repeats, angles_msg_dim_01, '-k')
xlabel('Milliseconds')
ylabel('$\|\sin\Theta\|_F^2$', 'interpreter', 'latex')
set(gca,'position',[0.1 0.1 0.85 0.43],'units','normalized')
axis([0, 5000, 0.1, 2])


%set(gca,'xtick',[])
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'--r');
h(2) = plot(NaN,NaN,'-k');
legend(h, 'red','blue','black');
legend(h, 'SGD', 'MSG')