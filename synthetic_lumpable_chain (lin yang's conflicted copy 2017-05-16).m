% lumpable states 3
% nodes 12

p0 = ...
[         0    0.1344    0.0463    0.0808    0.0463    0.1344    0.0463    0.0808    0.0808    0.1344    0.0808    0.1344
    0.1615         0    0.1615    0.0680    0.1615    0.0273    0.1615    0.0680    0.0680    0.0273    0.0680    0.0273
    0.0463    0.1344         0    0.0808    0.0463    0.1344    0.0463    0.0808    0.0808    0.1344    0.0808    0.1344
    0.1470    0.1030    0.1470         0    0.1470    0.1030    0.1470         0         0    0.1030         0    0.1030
    0.0463    0.1344    0.0463    0.0808         0    0.1344    0.0463    0.0808    0.0808    0.1344    0.0808    0.1344
    0.1615    0.0273    0.1615    0.0680    0.1615         0    0.1615    0.0680    0.0680    0.0273    0.0680    0.0273
    0.0463    0.1344    0.0463    0.0808    0.0463    0.1344         0    0.0808    0.0808    0.1344    0.0808    0.1344
    0.1470    0.1030    0.1470         0    0.1470    0.1030    0.1470         0         0    0.1030         0    0.1030
    0.1470    0.1030    0.1470         0    0.1470    0.1030    0.1470         0         0    0.1030         0    0.1030
    0.1615    0.0273    0.1615    0.0680    0.1615    0.0273    0.1615    0.0680    0.0680         0    0.0680    0.0273
    0.1470    0.1030    0.1470         0    0.1470    0.1030    0.1470         0         0    0.1030         0    0.1030
    0.1615    0.0273    0.1615    0.0680    0.1615    0.0273    0.1615    0.0680    0.0680    0.0273    0.0680         0 ]

%%%re generate
n1 = 60;
deg = 3;
p0 = createRandLumpGraph(n1, deg, 0.5)
n1 = length(p0)

figure
subplot(2,1,1)
p0 = diag(sum(p0, 2))\ p0;
Gp = digraph(p0);
zp = plot(Gp, 'k')
zp.EdgeCData = Gp.Edges.Weight;
colormap gray
colormap(flipud(colormap))
%colorbar
xticks('')
yticks('')
%pbaspect([2 1 1])
set(gca,'position',[0.01 0.51 0.88 0.48],'units','normalized')
axis([-2.4, 2.5, -1.7, 1.9])
%fig = gcf;
%fig.PaperPositionMode = 'auto'
%fig_pos = fig.PaperPosition;
%fig.PaperSize = [fig_pos(3) fig_pos(4)];
%c1 = colorbar
%print(fig,'lumpable_example_full','-dpdf')

%ylabel(c1, 'Transition Probability')

[us, ds] = eig(p0');
us = us / sum(us(:,1));
u0 = us(:,1);
pt = diag(u0) * p0;
pt = (pt + pt')/2;

[us1, ds1, vs1] = svd(pt)
ll = diag(ds1)
p1 = diag(u0) \ (ll(1) * us(:,1) * us(:,1)' + ll(2) * us(:,2) * us(:,2)' + ll(3) * us(:,3) * us(:,3)')

group = kmeans(diag(u0) \ us(:, 1:3), 3);
%zp.NodeLabel = group

pr = zeros(3,3);
for i=1:3
    for j=1:3
        pij = sum(p0((group == i), (group == j)), 1)
        pr(i, j) = pij(1);
    end
    pr(i, :) = pr(i, :) / sum(pr(i, :))
end

hold on
subplot(2,1,2)
Gpr = digraph(pr);
zpr = plot(Gpr, 'k')
zpr.EdgeCData = Gpr.Edges.Weight;
colormap gray
colormap(flipud(colormap))
%colorbar
xticks('')
yticks('')
%c1 = colorbar
%ylabel(c1, 'Transition Probability')
set(gca,'position',[0.01 0.01 0.88 0.48],'units','normalized')
axis([-1.4, 1.5, -1.2, 1.5])
zpr.NodeLabel=[{mat2str(find(group==1)')}, {mat2str(find(group==2)')}, {mat2str(find(group==3)')}]%{, find(group == 2), find(group == 3)}

h=colorbar;
caxis([0, 1])
set(h, 'position', [0.9, 0.1, 0.05, 0.8])


