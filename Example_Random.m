%Script  -- First Example, with random matrices in the manifold
clear all; close all; clc
rng(10)
n = 100;

%initialize the steady vector and its square root (entry-wise)
pi = rand(n,1); pi = pi/norm(pi,1);
pv = pi.^(1/2);

M = multinomialsymmetricfixedfactory(pv);

delta = 1e-5;
E = delta*rand(n);
Dv = diag(pv);

A = M.rand();

Ahat = diag(pv.^(-1))*A*diag(pv); %This is reversible, Ahat*e = e, pi'*Ahat = pi'
A0 = Ahat + E; %perturbation of Ahat

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
norm(diag(pi)*Xs - Xs'*diag(pi))
norm(Xs*ones(n,1) - ones(n,1))
norm(pi'*Xs - pi')

fprintf("||Xs - Ahat|| = %e\n", norm(Xs - A0,"fro"));
fprintf("Size initial pert. ||E|| = %e\n", norm(Ahat - A0,"fro"))

 infotable = struct2table(info);
% e = sqrt(infotable.cost);
% t = infotable.time;