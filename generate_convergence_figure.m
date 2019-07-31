figure
load('run_100_times.mat')
subplot(2,1,1)

lx =(1: 10^5);
init_exp = median(initial_all);
for rpt = 1:repeats
    semilogy((angle_all(rpt, 1: chain_length)), '--')
    hold on
end
semilogy((init_exp*exp(-abs((ll(k)-ll(k+1)) / ll(1)) * eta * lx / 12)), '-k')
axis([1, 7e4, 0.002, 10])
ylabel('$$\|\sin\Theta\|_F^2$$', 'Interpreter', 'latex')
%xlabel('Iterations')
set(gca,'xtick',[])
set(gca,'position',[0.1 0.55 0.85 0.43],'units','normalized')
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'--r');
h(2) = plot(NaN,NaN,'-k');
legend(h, 'red','blue','black');
legend(h, 'Algorithm', 'Thoeretic upper bound')

subplot(2,1,2)
load('run_100_times_diminishing.mat')
init_exp = median(initial_all);
for rpt = 1:repeats
    semilogy((angle_all(rpt, 1: chain_length)), '--')
    hold on
end
semilogy((init_exp*exp(-abs((ll(k)-ll(k+1)) / ll(1)) * eta * lx / 4)), '-k')
axis([1, 7e4, 0.002, 10])
set(gca,'position',[0.1 0.10 0.85 0.43],'units','normalized')
ylabel('$$\|\sin\Theta\|_F^2$$', 'Interpreter', 'latex')
xlabel('Iterations')

