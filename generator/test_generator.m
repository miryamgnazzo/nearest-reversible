%% Test Generator

clear; clc; close all;

classes = ["uniform","normal","sbm","multipleergodic"];
seed = 17;

sizes = repmat([50,100,200],[length(classes),1]);
number = repmat([5,5,5],[length(classes),1]);

[P,pi] = markov_generator(classes,sizes,number,seed);

% We save the test set to file, in this way we avoid generating it every
% time
save('test_set.mat','P','pi');

%% Visualize
figure(1)
heatmap(P{1,1})
set(gcf,'Color','white')
colorbar('off')
try
    export_fig("uniform_markov.pdf")
catch
end
figure(2)
heatmap(P{2,1})
set(gcf,'Color','white')
colorbar('off')
try
    export_fig("normal_markov.pdf")
catch
end
figure(3)
heatmap(P{3,1})
set(gcf,'Color','white')
colorbar('off')
try
    export_fig("sbm_markov.pdf")
catch
end
figure(4)
heatmap(P{4,1})
set(gcf,'Color','white')
colorbar('off')
try
    export_fig("mergodic_markov.pdf")
catch
end
