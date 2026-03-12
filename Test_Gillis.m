%Numerical test, employing the algorithm by Gillis and Van Dooren,
%available at https://gitlab.com/ngillis/TSDP 
clear all; close all; clc
rng(2)

n = 8;
d0 = rand(n,1);
d1 = rand(n-1,1); 
d2 = rand(n-2,1);

P = diag(d0) + diag(d1,1) + diag(d1,-1) + diag(d2,2) + diag(d2,-2);
P(1,n) = rand;
P(n,1) = P(1,n);
P = diag(sum(P,2))\P;

ps = rand(n,1); ps = abs(ps);
ps = ps/sum(ps);

A = zeros(n);

for i = 1:n
    for j = 1:n
        if i ~= j && P(i,j) > 0
            % Using Metropolis-Hastings formula
            mh_ratio = (ps(j) * P(j,i)) / (ps(i) * P(i,j));
            alpha = min(1, mh_ratio);
            A(i,j) = P(i,j) * alpha;
        end
    end
end
%Make the matrix stochastic
for i = 1:n
    A(i,i) = 1 - sum(A(i,[1:i-1, i+1:end]));
end

%Perturbation to make it full
Pert = 1e-5.*rand(n);

Atilde = A+Pert;
row_sums = sum(Atilde,2);
Atilde = diag(1./row_sums) * Atilde;

options.support = ones(n);
%Use Gurobi(default) in supportTSDP, NOT Matlab

% Apply Gillis and Van Dooren code
[Delta,~,~,~] = supportTSDP(Atilde,ps, options);
S = Atilde + Delta;

format short e
fprintf('Reversibility : %d \n', norm(diag(ps)*S - S'*diag(ps),'inf'));
fprintf('Stochastic : %d \n', norm(S*ones(n,1) - ones(n,1),'inf'));
fprintf('Stationary : %d\n', norm(ps'*S - ps','inf'));

[R,info] = riemannian_nearest_reversible(S,ps,'verbose',true,'RecurseErgodic',true);

norm(R-S,'fro')
norm(R-Atilde,'fro')
norm(R-A,'fro')

