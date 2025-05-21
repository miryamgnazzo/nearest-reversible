function [P,info] = riemannian_nearest_reversible(A,pi,varargin)
%RIEMANNIAN_NEAREST_REVERSIBLE: given a matrix A, computes the nearest reversible stochastic matrix P
%with the same stationary vector pi via Riemannian optimization
%The algorithm identifies transient states and possibly separates the problem into ergodic classes (if known)

%% Parsing of the inputs
p = inputParser;

addRequired(p,'A',@ismatrix)
addRequired(p,'pi',@iscolumn)
addOptional(p,'RecurseErgodic',false,@islogical)
addOptional(p,'verbose',false,@islogical)

parse(p,A,pi,varargin{:})

RecurseErgodic = p.Results.RecurseErgodic;
verbose = p.Results.verbose;

%% Start Computations
if (verbose) 
    fprintf('Verbose mode on\n');
end 

[n,m] = size(A);
% Check inputs
if n ~= m, error("riemannian_nearest_reversible:: is not square (%d,%d)",n,m), end
if norm(sum(P,2)-1,"inf") > 10*eps
    if any(sum(P,2) == 0)
        if (verbose); fprintf("The chain has 0 row-sums\n"); end
    else
        error('riemannian_nearest_reversible:: P is not Stochastic! norm(sum(P,2)-1,"inf") = %e > %e',norm(sum(P,2)-1,"inf"),10*eps);
    end
end
if any(pi < 0)
    if all(pi < 0)
        pi = -pi;
    else
        error('riemannian_nearest_reversible:: stationary distribution has changing signs')
    end
end

% Time variables:
time_presolve = 0;
time_solve = 0;

tic;
zeroind = find(abs(pi) < 100*eps);
if ~isempty(zeroind)
    if verbose, fprintf('P has transient states, reducing problem to recurrent subchain(s)'), end
    Q = A(setdiff(1:n,zeroind),setdiff(1:n,zeroind));
    qpi = pi(setdiff(1:n,zeroind),1);
else
    Q = A;
    qpi = pi;
end
time_presolve = time_presolve + toc;

if RecurseErgodic
    % We look for the different ergodic classes
    tic;
    if (verbose), fprintf('Looking for the existence of ergodic classes.\n'), end
    G = digraph(Q);
    [bins,binsize] = G.conncomp("OutputForm","vector","Type","strong");
    
    E = length(binsize); 
    if (verbose), fprintf('Found %d Ergodic Classes.\n',E), end
    QE = cell(E,1);
    piE = cell(E,1);
    REc = cell(E,1);

    for i = 1:E
        QE{i} = Q(bins == i,bins == i);
        piE{i} = qpi(bins == i,1);
    end
    time_presolve = time_presolve + toc;
    
    %% Solve phase
    % We run the solver on each ergodic class
    tic;
    for i = 1:E
        if length(QE{i}) > 1 
            if verbose, fprintf("Running Optimization on class %d of size %d x %d.\n",i,size(QE{i})),end
            % If the ergodic class is made by more than a single state

            %set up for manopt
            pv = piE{i}.^(1/2);
            M = multinomialsymmetricfixedfactory(pv);

            %clear problem; %If needed
            problem.M = M;
            problem.cost = @(X) 0.5*norm(diag(pv.^(-1))*X*diag(pv) - QE{i},'fro')^2;
            problem.egrad = @(X) (diag(piE{i}.^(-1))*X*diag(piE{i}) - diag(pv.^(-1))*QE{i}*diag(pv));
            problem.ehess = @(X, dX) diag(piE{i}.^(-1))*dX*diag(piE{i});
            
            options.verbosity = 0;

            [X1, xcost, ~, ~] = trustregions(problem, [], options);
            REc{i} = diag(pv.^(-1))*X1*diag(pv); %Solution of the problem - Reversible
        else
            if verbose, fprintf("Skipping solution, class is of size 1\n"), end
            % The ergodic class is an isolated state: we are reversible
            REc{i} = QE{i};
        end
    end
    time_solve = time_solve + toc;
    
    %% Assemble back the matrix - Phase 1
    % Assemble with respect to ergodic classes
    RE = zeros(size(Q));
    for i = 1:E
        RE(bins == i,bins == i) = REc{i};
    end
else
    % We don't care if there are any other ergodic classes and run the
    % Riemannian optimization algorithm on the whole Q
    if verbose, fprintf("Running Optimization on chain of size %d x %d.\n",size(Q)),end
   
    %set up for manopt
    pv = qpi.^(1/2);
    M = multinomialsymmetricfixedfactory(pv);
    problem.M = M;
    problem.cost = @(X) 0.5*norm(diag(pv.^(-1))*X*diag(pv) - Q,'fro')^2;
    problem.egrad = @(X) (diag(qpi.^(-1))*X*diag(qpi) - diag(pv.^(-1))*Q*diag(pv));
    problem.ehess = @(X, dX) diag(qpi.^(-1))*dX*diag(qpi);

    options.verbosity = 0;
    
    tic;
    [X1, xcost, ~, ~] = trustregions(problem, [], options);
    RE = diag(pv.^(-1))*X1*diag(pv); %Solution of the problem - Reversible
    time_solve = time_solve + toc;
end

%% Assemble back the matrix - Phase 2
% Assemble with respect to transient states
if ~isempty(zeroind)
    % We had transient states
    P = zeros(n,n);
    P(setdiff(1:n,zeroind),setdiff(1:n,zeroind)) = RE;
    P(zeroind,:) = A(zeroind,:);
else
    P = RE;
end

%% Output
if nargout > 1
    info.time_presolve = time_presolve;
    info.time_solve = time_solve;
    Dpi = spdiags(pi,0,n,n);
    info.reversibility = norm(Dpi*P - P'*Dpi,"inf");
    info.stochasticity = norm(sum(P,2)-1,"inf");
    info.stationarity = norm(pi'*P-pi',"inf");
    info.distance = norm(A-P,"fro");
    info.relative_distance = info.distance/norm(A,"fro");
end

end
