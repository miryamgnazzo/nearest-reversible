%% Test Generator

clear; clc; close all;

classes = ["uniform","normal","sbm"];
seed = 17;

sizes = repmat([50,100,200,300,500],[length(classes),1]);
number = repmat([20,20,20,20,20],[length(classes),1]);

[P,pi] = markov_generator(classes,sizes,number,seed);