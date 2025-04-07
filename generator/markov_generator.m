function [P,pi] = markov_generator(classes,sizes,number,seed)
%MARKOV_GENERATOR Generates number matrices of classes "class" and
%prescribed size.
%   INPUT:
%       class: string array denoting the type of Markov chain
%       size:  matrix of sizes size(i,:) contains the sizes for the ith
%       number: number of test matrices number(i,j) how many matrices of
%       class(i) and size size(i,j) to be generated
%       seed (optional): seed for the random generator if it is not given
%       is fixed to 1.
%   OUTPUT:
%       P cell array of size class x max(sum(number,2)) containing the
%       generated transition matrices
%       pi stationary distribution of the related matrix P

if exist("seed","var")
    rng(seed);
else
    rng(1);
end

P = cell(length(classes),max(sum(number,2)));
pi = cell(length(classes),max(sum(number,2)));

classind = 1;
for class = classes
    testm = 1;
    for size = sizes(classind,:)
        for num = number(classind,:)
            for j = 1:num
                [P{classind,testm},pi{classind,testm}] = genmarkov(class,size);
                testm = testm + 1;
            end
        end
    end
    classind = classind + 1;
end



end

function [P,pi] = genmarkov(class,Size)

switch lower(class)
    case "uniform"
        A = rand(Size,Size);
        D = sum(A,2);
        P = diag(D)\A;
    case "normal"
        A = 1 + randn(Size,Size);
        D = sum(A,2);
        P = diag(D)\A;
    case "sbm"
        isconnected = false;
        while ~isconnected
            % SBM parameters
            block_sizes = int32(get_random_block_sizes(Size));
            prob_matrix = generate_diagonally_dominant_matrix(length(block_sizes));
            % Convert MATLAB data to Python-compatible formats
            py_block_sizes = py.list(block_sizes);
            rows = size(prob_matrix, 1);
            py_prob_matrix = py.list();
            for i = 1:rows
                py_row = py.list(prob_matrix(i, :));
                py_prob_matrix.append(py_row);
            end
            % Call Python function
            adj_matrix_py = py.sbm_generator.generate_sbm(py_block_sizes, py_prob_matrix);
            % Convert the result back to MATLAB
            A = double(adj_matrix_py);
            D = sum(A,2);
            if all(D > 0)
                isconnected = true;
            end
        end
        P = diag(D)\A;
    case "multipleergodic"
        block_sizes = int32(get_random_block_sizes(Size));
        blocks = cell(length(block_sizes),1);
        for i=1:length(block_sizes)
            blocks{i} = rand(block_sizes(i),block_sizes(i));
        end
        A = blkdiag(blocks{:});
        D = sum(A,2);
        P = diag(D)\A;
    otherwise
        error("Unknown class %s",class)
end

[pi,~] = eigs(P',1,"largestabs");
pi = pi.*sign(pi);
pi = pi./sum(pi);
end