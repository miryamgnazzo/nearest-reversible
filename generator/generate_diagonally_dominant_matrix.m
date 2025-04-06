function P = generate_diagonally_dominant_matrix(n)
% Generate a random diagonally dominant matrix of size n x n
% All values are in the range [0, 1]

A = rand(n); % Random values in [0,1]

for i = 1:n
    % Sum of off-diagonal values in row i
    off_diag_sum = sum(A(i,:)) - A(i,i);

    % Ensure diagonal is at least as large as sum of off-diagonal values + epsilon
    epsilon = rand() * 0.1;  % Small positive value to ensure strict dominance
    A(i,i) = min(1, off_diag_sum + epsilon); % Clamp to 1 to stay in [0,1]
end

D = sum(A,2);
P = diag(D)\A;
end