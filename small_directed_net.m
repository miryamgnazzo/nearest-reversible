clear; clc; close all;

try
    fprintf("Manopt version: %d.%d.%d\n",manopt_version);
catch
    cd('manopt');
    addpath(genpath(pwd()));
    cd('..');
    fprintf("Manopt version: %d.%d.%d\n",manopt_version);
end

addpath('ia-southernwomen/');
[s,t] = importfile('ia-southernwomen.edges');
G = digraph(s,t,"omitselfloops");
figure(1)
subplot(1,2,1)
hG = plot(G,'NodeColor','k');

figure(1)
A = adjacency(G);
subplot(1,2,2);
spy(A);
xticklabels([])
yticklabels([])
xticks([])
yticks([])
set(gcf,'Color','white')

% Build Markov-Chain
D = diag(sum(A,2));
P = D\A;

% Compute Stationary Distribution
[pi,~] = eigs(P',1,'largestabs');
pi = pi.*sign(pi);
pi = pi./sum(pi);

% Find nearest reversible
[R,info] = riemannian_nearest_reversible(P,pi,'verbose',true,'RecurseErgodic',true);
disp(info)

% Visualize R
figure(2)
subplot(1,2,1)
heatmap(R)

subplot(1,2,2)
GR = digraph(R);
hGR = plot(GR,"XData",hG.XData,"YData",hG.YData,'NodeColor','k');
for i = 1:length(s)
    if GR.findedge(s(i),t(i)) && s(i) ~= t(i)
        highlight(hGR,s(i),t(i),'EdgeColor','k','LineStyle',':','LineWidth',0.5);
    end
end
xticklabels([])
yticklabels([])
xticks([])
yticks([])
xlabel(sprintf('nz = %d',nnz(R)))
axis square
set(gcf,'Color','white')

% Arc-weight
figure(3)
hGR = plot(GR,"XData",hG.XData,"YData",hG.YData,"LineWidth",GR.Edges.Weight*10,'NodeColor','k');
for i = 1:length(s)
    if GR.findedge(s(i),t(i)) && s(i) ~= t(i)
        highlight(hGR,s(i),t(i),'EdgeColor','k','LineStyle',':','LineWidth',P(s(i),t(i))*10);
    end
end
xticklabels([])
yticklabels([])
xticks([])
yticks([])
xlabel(sprintf('nz = %d',nnz(R)))
axis square
set(gcf,'Color','white')

%% QP approach
% Visualize R
R = getClosestSparse(P,pi);
figure(4)
subplot(1,2,1)
heatmap(R)

subplot(1,2,2)
GR = digraph(R);
hGR = plot(GR,"XData",hG.XData,"YData",hG.YData,'NodeColor','k');
for i = 1:length(s)
    if GR.findedge(s(i),t(i)) && s(i) ~= t(i)
        highlight(hGR,s(i),t(i),'EdgeColor','k','LineStyle',':','LineWidth',0.5);
    end
end
xticklabels([])
yticklabels([])
xticks([])
yticks([])
xlabel(sprintf('nz = %d',nnz(R)))
axis square
set(gcf,'Color','white')

% Arc-weight
figure(5)
hGR = plot(GR,"XData",hG.XData,"YData",hG.YData,"LineWidth",GR.Edges.Weight*10,'NodeColor','k');
for i = 1:length(s)
    if GR.findedge(s(i),t(i)) && s(i) ~= t(i)
        highlight(hGR,s(i),t(i),'EdgeColor','k','LineStyle',':','LineWidth',P(s(i),t(i))*10);
    end
end
xticklabels([])
yticklabels([])
xticks([])
yticks([])
xlabel(sprintf('nz = %d',nnz(R)))
axis square
set(gcf,'Color','white')

Dpi = diag(pi);
fprintf('QP reverse: %1.2e\n',norm(Dpi*R - R'*Dpi,"inf"));
fprintf('QP stationary: %1.2e\n',norm(pi'*R - pi',"inf"));
fprintf('QP distance: %1.2e\n',norm(R-P,"fro")/norm(P,"fro"));