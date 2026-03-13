%% Running the test on synthetic Markov chains

clear; clc; close all;

addpath('generator/');
try
    fprintf("Manopt version: %d.%d.%d\n",manopt_version);
catch
    cd('manopt');
    addpath(genpath(pwd()));
    cd('..');
    fprintf("Manopt version: %d.%d.%d\n",manopt_version);
end
addpath('/software/gurobi/gurobi1102/linux64/matlab');

%% Load test problems
load('test_set_smaller.mat');
fprintf("Number of classes: %d\n",size(P,1));
fprintf("Number of test problems: %d x %d = %d\n",size(P,1),size(P,2),numel(P));

total_test_number = numel(P);
number_of_algorithms = 3;

solve_time     = NaN(total_test_number,number_of_algorithms);
distance_norm  = NaN(total_test_number,number_of_algorithms);
reversibility  = NaN(total_test_number,number_of_algorithms);
stationary     = NaN(total_test_number,number_of_algorithms);
stochastic     = NaN(total_test_number,number_of_algorithms);
actual_count = 0;

for i = 1:total_test_number
    fprintf("Test number: %d\n",i)

    Pmat = P{i};
    pivec = pi{i};
    Dpi = diag(pivec);

    % Check if Pmat has ergodic classes
    G = digraph(Pmat);
    [bins,binsize] = G.conncomp("OutputForm","vector","Type","strong");
    hasergodic = length(binsize) > 1;
    
    if size(Pmat,1) > 0
        
        actual_count = actual_count + 1;

        %% Quadratic programming solver
        solver_number = 1;
        try
            fprintf("Solving problem %d with MATLAB QP solver\n",i);
            [R,info] = qp_nearest_reversible(Pmat,pivec,...
                'RecurseErgodic',hasergodic,'verbose',false,'solver','quadprog');
            solve_time(i,solver_number) = info.time_presolve + info.time_solve;
            distance_norm(i,solver_number) = info.relative_distance;
            reversibility(i,solver_number) = info.reversibility;
            stationary(i,solver_number) = info.stationarity;
            stochastic(i,solver_number) = norm(1-sum(R,2),"inf");
        catch
            fprintf("Solver %d failed on problem %d :-(\n",solver_number,i);
        end
        

        %% Riemannian solver
        solver_number = 2;
        try
            fprintf("Solving problem %d with Riemann solver\n",i);
            [R,info] = riemannian_nearest_reversible(Pmat,pivec,...
                'RecurseErgodic',true,'verbose',false);
            solve_time(i,solver_number) = info.time_presolve + info.time_solve;
            distance_norm(i,solver_number) = info.relative_distance;
            reversibility(i,solver_number) = info.reversibility;
            stationary(i,solver_number) = info.stationarity;
            stochastic(i,solver_number) = norm(1-sum(R,2),"inf");
        catch
            fprintf("Solver %d failed on problem %d :-(\n",solver_number,i);
        end

        %% Gurobi solver (Automatic)
        solver_number = 3;
        try
            fprintf("Solving problem %d with GUROBI QP solver\n",i);
            [R,info] = qp_nearest_reversible(Pmat,pivec,...
                'RecurseErgodic',hasergodic,'verbose',false,'solver','gurobi');
            solve_time(i,solver_number) = info.time_presolve + info.time_solve;
            distance_norm(i,solver_number) = info.relative_distance;
            reversibility(i,solver_number) = info.reversibility;
            stationary(i,solver_number) = info.stationarity;
            stochastic(i,solver_number) = norm(1-sum(R,2),"inf");
        catch
            fprintf("Solver %d failed on problem %d :-(\n",solver_number,i);
        end
               
        %% Gurobi solver (Primal Simplex)
        % solver_number = 4;
        % try
        %    [R,solve_time(i,solver_number)] = getClosestSparse_gurobi(Pmat,pivec,0);
        %    distance_norm(i,solver_number) = norm(Pmat-R,"fro")/norm(Pmat,"fro");
        %    reversibility(i,solver_number) = norm(Dpi*Pmat - Pmat'*Dpi,"inf");
        %    stationary(i,solver_number) = norm(pivec'*Pmat-pivec',"inf");
        % catch
        %    fprintf("Solver %d failed on problem %d :-(",solver_number,i);
        % end

    else
        fprintf("Skipped test: the matrix is empty!")
    end

    %% Save results to file
    save("test_runner_results_new.mat","stationary","reversibility","solve_time","distance_norm","stochastic")


end
