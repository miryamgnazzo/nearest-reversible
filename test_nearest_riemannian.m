%% Functionality Test for the Riemannian Nearest Reversible Implementation

clear; clc; close all;
addpath('generator/');

%% First we test a case without asking for Ergodic classes
load('test_set.mat');

[R,info] = riemannian_nearest_reversible(P{4,6},pi{4,6},'RecurseErgodic',...
    false,'verbose',true);

figure(1)
subplot(1,2,1); spy(P{4,6});
subplot(1,2,2); spy(R);


%% But this test case has two, so we look for classes
[R,info] = riemannian_nearest_reversible(P{4,6},pi{4,6},'RecurseErgodic',...
    true,'verbose',true);

figure(2)
subplot(1,2,1); spy(P{4,6});
subplot(1,2,2); spy(R);