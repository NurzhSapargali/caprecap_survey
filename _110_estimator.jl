module Estimator

import Distributions: Beta

using NLopt
using StatsBase

export loglh, fit_model


function monte_carlo(prob_matrix::Matrix{Float64},
                     complement_prob_matrix::Matrix{Float64},
                     x_i::Vector{Bool})
    X = repeat(x_i, 1, size(prob_matrix)[2]);
    G_X = (prob_matrix .^ X) .* (complement_prob_matrix .^ (1 .- X));
    integrands = prod(G_X, dims=1);
    return mean(integrands);
end

function loglh(alpha::Float64,
               N_u::Float64,
               S::Vector{Vector{Int64}},
               O::Set{Int64},
               n::Vector{Int64},
               draws::Int64)
    N_o = length(O);
    T = length(S);
    beta = alpha * (N_u + N_o - 1.0);
    points = rand(Beta(alpha, beta), draws);
    P = n * transpose(points);
    comp_P = 1.0 .- P;
    sum_term = 0.0;
    for i in O
        x_i = [i in s for s in S];
        I = monte_carlo(P, comp_P, x_i);
        if I < 0
            I = 5e-200;
        end
        sum_term += log(I);
    end
    truncation = 1.0 - monte_carlo(P, comp_P, zeros(Bool, T));
    if truncation < 0
        truncation = 5e-200;
    end
    lh = -N_o * log(truncation) + sum_term;
    println("alpha = $alpha, N_u = $N_u, lh = $lh");
    return lh;
end

function fit_model(S::Vector{Vector{Int64}},
                   O::Set{Int64},
                   n::Vector{Int64},
                   draws::Int64)
    LL(x, grad) = -loglh(x[1], x[2], S, O, n, draws);
    opt = Opt(:LN_SBPLX, 2);
    lower = [0.01, 0];
    upper = [10000, Inf];
    opt.upper_bounds = upper
    opt.lower_bounds = lower;
    opt.min_objective = LL;
    (minf, minx, ret) = NLopt.optimize(opt, [5.0, length(O)]);
    return (minf, minx, ret);
end

# function simulated_annealing(accepts, sigma_alpha, sigma_Nu,
#                              alpha_0, N_u_0, S, O, n, mc_draws)
#     alphas = Vector{Float64}();
#     N_us = Vector{Float64}();
#     evals = Vector{Float64}();
#     temps = Vector{Float64}();
#     run = 1.0;
#     alpha_old = alpha_0;
#     N_u_old = N_u_0;
#     while length(evals) < accepts
#         T = 1.0 / log(1.0 + run);
#         lh_old = loglh(alpha_old, N_u_old, S, O, n, mc_draws);
#         alpha_new = alpha_old + rand(Normal(0.0, sigma_alpha));
#         N_u_new = N_u_old + rand(Normal(0.0, sigma_Nu));
#         if (N_u_new < 0)
#             N_u_new = 0.0
#         end
#         if (alpha_new <= 0.0)
#             alpha_new = 1e-3
#         end
#         lh_new = loglh(alpha_new, N_u_new, S, O, n, mc_draws);
#         del_lh = lh_new - lh_old;
#         rho = minimum([exp(del_lh / T), 1.0]);
#         rho = maximum([rho, 0.0]);
#         accept = rand(Bernoulli(rho));
#         if accept
#             alpha_old = alpha_new;
#             N_u_old = N_u_new;
#             lh_new = loglh(alpha_new, N_u_new, S, O, n, mc_draws);
#             push!(alphas, alpha_old);
#             push!(N_us, N_u_old);
#             push!(evals, lh_new);
#             push!(temps, T);
#             println("alpha = $alpha_old, N_u = $N_u_old, lh = $lh_new");
#         end
#         run += 1;
#     end
#     df = DataFrame(alpha = alphas, N_u = N_us, lh = evals, temp = temps);
#     return df;
# end

# function spsa(iterations, alpha_0, N_u_0, a,
#               c, gamma, S, O,
#               n, mc_draws)
#     counter = 1.0;
#     alpha_k = alpha_0;
#     N_u_k = N_u_0;
#     while counter < iterations
#         a_k = a / counter;
#         c_k = c / (counter ^ gamma);
#         delta_k = [-1.0, -1.0] .^ rand(Bernoulli(0.5));
#         pert_plus = [alpha_k, N_u_k] + c_k * delta_k;
#         if pert_plus[1] <= 0.0
#             pert_plus[1] = 1e-3; 
#         end
#         if pert_plus[2] < 0.0
#             pert_plus[2] = 0.0; 
#         end
#         pert_minus = [alpha_k, N_u_k] - c_k * delta_k;
#         if pert_minus[1] <= 0.0
#             pert_minus[1] = 1e-3; 
#         end
#         if pert_minus[2] < 0.0
#             pert_minus[2] = 0.0; 
#         end
#         lh_plus = -loglh(pert_plus[1], pert_plus[2], S, O, n, mc_draws);
#         lh_minus = -loglh(pert_minus[1], pert_minus[2], S, O, n, mc_draws);
#         grad = (lh_plus - lh_minus) / (2.0 * c_k) .* (1.0 ./ delta_k);
#         alpha_k = alpha_k - a_k * grad[1];
#         if alpha_k <= 0.0
#             alpha_k = 1e-3;
#         end
#         N_u_k = N_u_k - a_k * grad[2];
#         if N_u_k < 0.0
#             N_u_k = 0.0;
#         end
#     end
#     return [alpha_k, N_u_k];
# end

end
