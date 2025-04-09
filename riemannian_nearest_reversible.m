function [R,info] = riemannian_nearest_reversible(P,pi,varargin)
%RIEMANNIAN_NEAREST_REVERSIBLE Computes the nearest reversible Markov chain
%with the same stationary vector pi via Riemannian optimization
%   Detailed explanation goes here

%% Parsing of the inputs
p = inputParser;

addRequired(p,'P',@ismatrix)
addRequired(p,'pi',@iscolumn)
addOptional(p,'RecurseErgodic',false,@islogical)
addOptional(p,'verbose',false,@islogical)

parse(p,P,pi,varargin{:})

RecurseErgodic = p.Results.RecurseErgodic;
verbose = p.Results.verbose;

%% Start Computations
if (verbose) 
    fprintf('Verbose mode on\n');
end 

[n,m] = size(P);
% Check inputs
if n ~= m, error("riemannian_nearest_reversible:: is not square (%d,%d)",n,m), end
if norm(sum(P,2)-1,"inf") > 10*eps
    error('riemannian_nearest_reversible:: P is not Stochastic! norm(sum(P,2)-1,"inf") = %e > %e',norm(sum(P,2)-1,"inf"),10*eps);
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
zeroind = find(abs(pi) < 10*eps);
if ~isempty(zeroind)
    if verbose, fprintf('P has transient states, reducing problem to recurrent subchain(s)'), end
    Q = P(setdiff(1:n,zeroind),setdiff(1:n,zeroind));
else
    Q = P;
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
        piE{i} = pi(bins == i,1);
    end
    time_presolve = time_presolve + toc;
    
    %% Solve phase
    % We run the solver on each ergodic class
    tic;
    for i = 1:E
        if length(QE{i}) > 1 
            if verbose, fprintf("Running Optimization on class %d of size %d x %d.\n",i,size(QE{i})),end
            % If the ergodic class is made by more than a single state
            REc{i} = ones(size(QE{i})); %TODO: Sostituire qui l'ottimizzazione Riemanniana!
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
    tic;
    RE = ones(size(Q)); %TODO: Sostituire qui l'ottimizzazione Riemanniana!
    time_solve = time_solve + toc;
end

%% Assemble back the matrix - Phase 2
% Assemble with respect to transient states
if ~isempty(zeroind)
    % We had transient states
    R = zeros(n,n);
    R(setdiff(1:n,zeroind),setdiff(1:n,zeroind)) = RE;
    R(zeroind,zeroind) = P(zeroind,zeroind);
else
    R = RE;
end

%% Output
if nargout > 1
    info.time_presolve = time_presolve;
    info.time_solve = time_solve;
    Dpi = spdiags(pi,0,n,n);
    info.reversibility = norm(Dpi*R - R'*Dpi,"inf");
    info.stochasticity = norm(sum(R,2)-1,"inf");
    info.stationarity = norm(pi'*R-pi',"inf");
    info.distance = norm(P-R,"fro");
    info.relative_distance = info.distance/norm(P,"fro");
end

end