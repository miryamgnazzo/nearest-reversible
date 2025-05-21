%% Evaluation of memory consumption for the different algorithms

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

%% Generate uniform test problems
%[P,pi] = markov_generator("uniform",[10 100 200 300 500],[1 1 1 1 1]);
load('matrix-memory.mat')

peakMem           = zeros(5,5,2);
TotalMemAllocated = zeros(5,5,2);



%% Start Profiling for the ClosestSparse routine
index = 1;
for i = 1:5
    for j=1:5
        profile on -memory  % Start profiling with memory tracking
        U = getClosestSparse(P{index}, pi{index});  % Call the function
        profile off  % Stop profiling
        p = profile('info');
        [peakMem(i,j,1),TotalMemAllocated(i,j,1)] = getmemory(p,'getClosestSparse');
        index = index + 1;
        clear p
    end
end

%% Start Profiling for the Riemannian routine
index = 1;
for i = 1:5
    for j=1:5
        profile on -memory  % Start profiling with memory tracking
        U = riemannian_nearest_reversible(P{index}, pi{index},...
            'RecurseErgodic',false,'verbose',false);  % Call the function
        profile off  % Stop profiling
        p = profile('info');
        [peakMem(i,j,2),TotalMemAllocated(i,j,2)] = getmemory(p,'riemannian_nearest_reversible');
        index = index + 1;
        clear p
    end
end

%% Save the output results to file
save("memory_usage.mat","P","pi","peakMem","TotalMemAllocated");

%% Auxiliary functions

function [peakMem,TotalMemAllocated] = getmemory(p,name)
% Initialize peak memory variable
peakMem = NaN;
TotalMemAllocated = NaN;
% Loop through function table to find 'name'
for k = 1:length(p.FunctionTable)
    funcInfo = p.FunctionTable(k);
    if contains(funcInfo.FunctionName, name)
        peakMem = funcInfo.PeakMem;
        TotalMemAllocated = funcInfo.TotalMemAllocated;
        fprintf('Peak memory used by %s: %.2f MB\n',name, peakMem / (1024^2));
        fprintf('TotalMemAllocated by %s: %.2f MB\n',name, TotalMemAllocated/(1024^2));
        break;
    end
end

% Optional: Handle case if function is not found
if isnan(peakMem) || isnan(TotalMemAllocated)
    warning('Function %s not found in profile data.',namedPattern);
end

end
