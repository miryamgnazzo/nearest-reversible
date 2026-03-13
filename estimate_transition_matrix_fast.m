function A_noisy = estimate_transition_matrix_fast(chain)
    % Estimate the transition matrix of a Markov chain
    K = max(chain);
    n = length(chain);
    
    % Extract 'from' states and 'to' states
    from_states = chain(1:end-1);
    to_states   = chain(2:end);
    
    % Count each transition
    C = full(sparse(from_states, to_states, 1, K, K));
    
    % Normalize to have a stochastic matrix
    row_sums = sum(C, 2);
    A_noisy = diag(1./row_sums)*C;
end