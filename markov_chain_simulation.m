function chain = markov_chain_simulation(A, nsteps, initial_state)
% MARKOV_CHAIN_SIMULATION: Simulates a discrete Markov chain.

% INPUT:
%   A:         K x K transition matrix.
%   nsteps:    Number of transition steps (total steps).
%   initial_state: The starting state (integer from 1 to K).
%
% OUTPUT:
%   chain:          Vector of length n+1 containing the sequence of states.

    K = size(A, 1);
    chain = zeros(1, nsteps + 1); 

  % Compute Cumulative Distribution Function (CDF)
    cdf_tot = cumsum(A, 2); 
    
    
  % Set the initial state
    chain(1) = initial_state;
    
  % Run the simulation loop
    for t = 1 : nsteps
        current_state = chain(t);
        cdf = cdf_tot(current_state,:);
        u = rand();
        next_state = find(cdf >= u, 1, 'first');
        chain(t + 1) = next_state;
    end
end