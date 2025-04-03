%Script -- Example for reducible matrices on a small problem
% Alcuni commenti: mi sembra che la soluzione perda di precisione e questo
% è dovuto alla lenta convergenza di SinkHorn (se ci sono entrate nulle)
% (ho dato anche un'occhiata qui:
% Philip A. Knight, "The Sinkhorn–Knopp Algorithm: Convergence and 
% Applications" in SIAM Journal on Matrix Analysis and Applications 30(1), 
% 261-275, 2008. )
% Per il momento ho aumentato le iterazioni di SinkHorn, ma mi sembra che
% il risultato di Neilsen-Weber sia comunque migliore
clear all; close all; clc
rng(1)
n = 10;
%initialize two different vectors
p1 = rand(n,1); p1 = p1/norm(p1,1);
p2 = rand(n,1); p2 = p2/norm(p2,1);

A1 = genrand(p1); A2 =  genrand(p2);
A = blkdiag(A1, A2);

a = 0.2;
pp = a*[p1; zeros(n,1)] + (1-a)*[zeros(n,1); p2];
pv = pp.^(1/2);

M = multinomialsymmetricfixedfactory(pv);

A0 = A;

%Initialization of the problem
problem.M = M;
problem.cost = @(X) cnormsqfro(diag(pv.^(-1))*X*diag(pv) - A0);
problem = manoptAD(problem);

% checkgradient(problem)

% options.tolgradnorm = 1e-10;
options.verbosity = 0;

%[X1, xcost, info, options] = conjugategradient(problem, [], options);
[X1, xcost, info, options] = trustregions(problem, [], options);

Xs = diag(pv.^(-1))*X1*diag(pv); %Solution of the problem

%Additional checks on the Reversibility
norm(diag(pp)*Xs - Xs'*diag(pp))
norm(Xs*ones(2*n,1) - ones(2*n,1))
norm(pp'*Xs - pp')

fprintf("||Xs - Ahat|| = %e\n", norm(Xs - A0,"fro"));

 infotable = struct2table(info);
% e = sqrt(infotable.cost);
% t = infotable.time;

