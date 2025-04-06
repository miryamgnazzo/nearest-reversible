function block_sizes = get_random_block_sizes(n)
    % Randomly choose the number of blocks between 2 and floor(n/2)
    k = randi([2, max(2, floor(n/4))]);
    
    % Randomly partition n into k parts (each >= 1)
    % Idea: choose k-1 unique cut points from 1 to n-1
    cuts = sort(randperm(n-1, k-1));
    cuts = [0, cuts, n];
    
    block_sizes = diff(cuts); % This gives k values that sum to n
end