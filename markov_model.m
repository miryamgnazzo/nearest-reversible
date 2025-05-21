%% Euler-Maruyama Method for SDE with 2pi-Periodic Potential and Markov Model

clear; clc; close all;

dt = 0.0001;    % Time step
N  = 6e7;      % Number of steps
x0 = 0;        % Initial condition
sigma = 1;     % Noise intensity
bin_count = 30; % Number of bins for Markov model
chunk_size = 400; % Subtrajectory length

Time = N*dt;  % Real time
x = zeros(N,1);
x(1) = x0;

fprintf("Running for %1.2e steps of size %1.2e to %1.2e.\n",N,dt,Time);

for i = 2:N
    dB = randn;  % Brownian motion increment
    x(i) = x(i-1) - potential_gradient(x(i-1)) * dt + sigma * sqrt(dt) * dB;
    x(i) = mod(x(i), 2*pi); % Apply periodic boundary condition
end

%% Visualize

time = linspace(0, Time, N);
figure(1)
plot(downsample(time(1:N/4),1000), downsample(x(1:N/4),1000));
xlabel('Time');
ylabel('$X_t$  mod  $2\pi$','Interpreter','LaTeX');
title('Euler-Maruyama Simulation of SDE with Periodic Potential');
bins = linspace(0, 2*pi, bin_count+1);
ylim([0,2*pi])
hline(bins,'k--')


%% Compute transition matrix
[C, T] = compute_markov_matrix(x, bin_count, chunk_size);
figure(1)
subplot(1,3,1)
heatmap(C)
title('Transition Count Matrix C:');
subplot(1,3,2)
heatmap(T)
title('Transition Probability Matrix T:');
subplot(1,3,3)
spy(T)

%% Check if T is reversible
[p,l] = eigs(T',1,"largestabs");
p = p.*sign(p); 
p = p./sum(p); % Normalize Steady State
p(abs(p) < 10*eps) = 0;
figure(2)
semilogy(1:length(T),p,'o')
title("Stationary distribution")

Dp = diag(p);
fprintf("||D T - T'*D||_inf = %e\n",norm(Dp,"inf"))

ev = eig(T);
ev = sort(ev,"descend");
figure(3)
plot(1:length(T),ev,'x')
xlabel("$k$",'Interpreter','latex')
ylabel("$\lambda_k$",'Interpreter','latex')


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