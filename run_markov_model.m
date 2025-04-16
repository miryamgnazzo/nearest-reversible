%% Euler-Maruyama Method for SDE with 2pi-Periodic Potential and Markov Model

clear; clc; close all;

%% Add Manopt
try
    fprintf("Manopt version: %d.%d.%d\n",manopt_version);
catch
    cd('manopt');
    addpath(genpath(pwd()));
    cd('..');
    fprintf("Manopt version: %d.%d.%d\n",manopt_version);
end

%% Status of the random number generator
% Set this to get reproducible results
rng(42);

%% Discretization parameters
dt = 0.001;       % Time step
N  = 5e8;         % Number of steps
x0 = 0;           % Initial condition
sigma = 1;        % Noise intensity
bin_count = 30;   % Number of bins for Markov model
chunk_size = 400; % Subtrajectory length

Time = N*dt;      % Real time
x = zeros(N,1);
x(1) = x0;

fprintf("Running for %1.2e steps of size %1.2e to %1.2e.\n",N,dt,Time);

%% Experimental parameters
nrepetition = 50;
time = NaN(2,nrepetition);
reversibility = NaN(2,nrepetition);
stationarity = NaN(2,nrepetition);
perturbation = NaN(2,nrepetition);
stochasticity = NaN(2,nrepetition);

for repetition = 1:nrepetition
    fprintf("Repetition n° %d\n",repetition)

    % Simulation with Euler-Maruyama method
    for i = 2:N
        dB = randn;  % Brownian motion increment
        x(i) = x(i-1) - potential_gradient(x(i-1)) * dt + sigma * sqrt(dt) * dB;
        x(i) = mod(x(i), 2*pi); % Apply periodic boundary condition
    end

    % Build Markov model
    [C, T] = compute_markov_matrix(x, bin_count, chunk_size);

    % Compute the steady-state vector
    [p,l] = eigs(T',1,"largestabs");
    p = p.*sign(p);
    p = p./sum(p); % Normalize Steady State
    p(abs(p) < 10*eps) = 0;
    Dp = diag(p);

    % Make it reversible the QP-way
    try
        tic;
        U = getClosestSparse(T, p);
        time(1,repetition) = toc;
        stationarity(1,repetition) = norm(p'*U - p',"inf");
        reversibility(1,repetition) = norm(Dp*U - U'*Dp,"inf");
        perturbation(1,repetition) = norm(U-T,"fro")/norm(T,"fro");
        stochasticity(1,repetition) = norm(sum(U,2)-1,"inf");
    catch
        fprintf("QP failed on repetition %d\n",repetition);
    end
    % Make it reversible the Riemannian-way
    try
        [R,info] = riemannian_nearest_reversible(T,p,"RecurseErgodic",true,"verbose",false);
        time(2,repetition) = info.time_presolve+info.time_solve;
        stationarity(2,repetition) = info.stationarity;
        reversibility(2,repetition) = info.reversibility;
        perturbation(2,repetition) = info.relative_distance;
        stochasticity(2,repetition) = info.stochasticity;
    catch
        fprintf("Riemann failed on repetition %d\n",repetition);
    end
end

%% Visualization as boxplots


labels = {'QP','Riemannian-Solver'};

figure(1)
boxplot(time','Labels',labels)
title('Time (s)')
matlab2tikz('filename','butan_time.tex','width','0.25\columnwidth');

figure(2)
boxplot(stationarity','Labels',labels)
hline(eps,'k--')
set(gca,'YScale','log')
title('Error on stationary distribution')
matlab2tikz('filename','butan_stationarity.tex','width','0.25\columnwidth');


figure(3)
boxplot(reversibility','Labels',labels)
hline(eps,'k--')
set(gca,'YScale','log')
title('Reversibility condition')
matlab2tikz('filename','butan_reversibility.tex','width','0.25\columnwidth');

figure(4)
boxplot(perturbation','Labels',labels)
set(gca,'YScale','log')
title('Relative distance')
matlab2tikz('filename','butan_perturbation.tex','width','0.25\columnwidth');

figure(5)
boxplot(perturbation','Labels',labels)
hline(eps,'k--')
set(gca,'YScale','log')
title('Relative distance')
matlab2tikz('filename','butan_stochasticity.tex','width','0.25\columnwidth');



%% Auxiliary functions
function gradV = potential_gradient(x)
% Define the gradient of the given periodic potential function V_B(x)
a = 2.0567;
b = -4.0567;
c = 0.3133;
d = 6.4267;

gradV = -b * sin(x) - 2 * c * cos(x) .* sin(x) - 3 * d * cos(x).^2 .* sin(x);
end

function [C, T] = compute_markov_matrix(x, bin_count, chunk_size)
% Compute the transition matrix from trajectory data
bins = linspace(0, 2*pi, bin_count+1); % Discretize state space
bin_indices = discretize(x, bins); % Assign states to bins

C = zeros(bin_count, bin_count); % Transition count matrix
num_chunks = floor(length(x) / chunk_size);

for k = 1:num_chunks-1
    i = bin_indices(k*chunk_size); % Start state
    j = bin_indices(k*chunk_size + 1); % Next state
    if ~isnan(i) && ~isnan(j)
        C(i,j) = C(i,j) + 1;
    end
end

% Normalize to get transition probability matrix
T = C ./ sum(C, 2);
T(isnan(T)) = 0; % Handle divisions by zero
end