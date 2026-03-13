%Test for Noisy experiment
clear all; close all; clc
rng(12)

n = 10;
P = rand(n,n);
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

%A is a reversible matrix

%Estimation of the Markov Chain from Nsteps sample
Nsteps = 100;
max_ex = 7;

obj_noisy = zeros(1,max_ex);
reconst_noisy = zeros(1,max_ex);

obj_exact = zeros(1,max_ex);
reconst_exact = zeros(1,max_ex);

norms = zeros(1,max_ex);
lower = zeros(1,max_ex);

for i = 1:max_ex

    X_0 = 1; %Start from the State 1
    
    chain = markov_chain_simulation(A, Nsteps, X_0);
    A_noisy = estimate_transition_matrix_fast(chain);
    
    [ps_noisy,~] = eigs(A_noisy',1,'largestabs');
    ps_noisy = ps_noisy.*sign(ps_noisy);
    ps_noisy = ps_noisy./sum(ps_noisy);
    
    
    [R_noisy,info_n] = riemannian_nearest_reversible(A_noisy,ps_noisy,'verbose',true,'RecurseErgodic',true);
    obj_noisy(i) = norm(R_noisy - A_noisy,'fro');
    reconst_noisy(i) = norm(R_noisy - A,'fro');
    
    [R_exact,info] = riemannian_nearest_reversible(A_noisy,ps,'verbose',true,'RecurseErgodic',true);
    obj_exact(i) = norm(R_exact - A_noisy,'fro');
    reconst_exact(i) = norm(R_exact - A,'fro');

    norms(i) = norm(R_exact-R_noisy,'inf');
    lower(i) = norm(ps-ps_noisy,'inf')/norm(ps_noisy,'inf');
    lower(i) = lower(i)/(norm(inv(eye(n)-R_exact + ones(n,1)*ps'),'inf'));

    Nsteps = Nsteps*10;
end

figure(1)
loglog(10.^[2:8], obj_noisy,'-*')
hold on
loglog(10.^[2:8], obj_exact,'-*')
loglog(10.^[2:8], reconst_noisy,'-*')
loglog(10.^[2:8], reconst_exact,'-*')

legend('Obj noisy', 'Obj Exact', 'Recon Noisy', 'Recon Exact')

figure(2)
loglog(10.^[2:8], norms,'-*')
hold on
loglog(10.^[2:8], lower,'-*')

legend('Norms', 'Lower Bound')