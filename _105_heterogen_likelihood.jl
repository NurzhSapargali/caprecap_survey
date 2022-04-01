using StatsBase
using Distributions
using Plots
 
n = [5; 10; 30; 5; 10; 20; 15; 5; 25; 60]; # Vector of sample sizes at each draw
N = 100; # Population size
ALPHA = 1.0; # Hyperparameter of prior


function loglike(Nu, O, s)
    product_term = 1;
    sum_term = 0;
    No = length(O);
    for t in 1:length(s)
        product_term *= 1 - length(t) / (Nu + No);
        for i in O
            x = 0;
            if i in s[t]
                x = 1;
            end
            sum_term += log(x * length(t) / (Nu + No) + (1 - x) * (1 - length(t) / (Nu + No)));
        end
    end
    return -No * log(1 - product_term) + sum_term;
end
w = ProbabilityWeights(rand(Beta(ALPHA, ALPHA * (N - 1)), N));
s = [sample(1:N, w, n[i], replace=false) for i in 1:length(n)];
O = Set();
for t in s
    O = union(O, Set(t));
end
y = [loglike(i, O, s) for i in 1:N];
plot(collect(1:N) , y)
print(argmax(y))
for experiment in 1:5
    w = ProbabilityWeights(rand(Beta(ALPHA, ALPHA * (N - 1)), N));
    s = [sample(1:N, w, n[i], replace=false) for i in 1:length(n)];
    O = Set();
    for t in s
        O = union(O, Set(t));
    end
    y = [loglike(i, O, s) for i in 1:N];
    print(argmax(y))
    plot!(collect(1:N), y)
end
