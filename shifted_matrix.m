%Test for matrices with zero entries - shifted version
clear all; close all; clc

classes = "multipleergodic";
seed = 11;
[P,pi] = markov_generator(classes, 50,1,seed);

P = P{1};
pi = pi{1};

a = 10.^(-[1:10]);
a = [a, 0.5, 0.7, 1];


s = length(a);
n = size(P,1);

e = ones(n,1);
negative = false(s,1);
distance = zeros(s,1);
rever = zeros(s,1);
stoch = zeros(s,1);
stat = zeros(s,1);
minval = zeros(s,1);

for i = 1:s

    Pa = a(i)*P + (1-a(i))*(e*pi');

    [R,info] = riemannian_nearest_reversible(Pa,pi,'RecurseErgodic',...
    false,'verbose',true);

    Q = (1/a(i))*(R - (1-a(i))*(e*pi'));

    negative(i) = any(Q(:) < 0); %check if there are any negative entries

    distance(i) = norm(P-Q,'inf');
    rever(i) = norm(diag(pi)*Q - Q'*diag(pi),"inf");
    stoch(i) = norm(sum(Q,2)-1,"inf");
    stat(i) = norm(pi'*Q-pi',"inf");
    minval(i) = min(min(Q));
end

T = table(a', negative, distance, rever, stoch, stat, minval);
disp(T)

