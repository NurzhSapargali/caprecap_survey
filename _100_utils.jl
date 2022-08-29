module Utils

import Distributions: Beta, Normal, Bernoulli
using StatsBase
using DataFrames

export sampford_sample, lincoln, schnabel,
       chao, chao_corrected, jackknife,
       loglh, simulated_annealing, write_row
       
function write_row(filename, row)
    open(filename, "a") do io
        for (i,j) in enumerate(row)
            if i != length(row)
                print(io, j, ",");
            else
                print(io, j, "\n");
            end
        end
    end
end

function sampford_sample(p, n)
    pr = n * p;
    pr = Dict(enumerate(pr));
    s = Set();
    while length(s) != n
        key = sample(1:length(p), Weights(p));
        val = pop!(pr, key);
        new_pr = collect(values(pr));
        new_pr = new_pr ./ (1.0 .- new_pr);
        new_pr = new_pr ./ sum(new_pr);
        s = Set(sample(collect(keys(pr)), Weights(new_pr), n));
        pr[key] = val;
    end
    println("Generated")
    return collect(s);
end

function monte_carlo(prob_matrix, complement_prob_matrix, x_i)
    X = repeat(x_i, 1, size(prob_matrix)[2]);
    G_X = (prob_matrix .^ X) .* (complement_prob_matrix .^ (1 .- X))
    integrands = prod(G_X, dims=1);
    return mean(integrands)
end

function lincoln(S, n)
    r = intersect(Set(S[1]), Set(S[2]));
    return prod(n) / length(r);
end

function schnabel(S, n)
    m = 0;
    N = 0;
    pool = Set();
    denom = [];
    num = [];
    R = Set();
    for t in 2:length(S)
        pool = union(pool, S[t-1]);
        u = n[t-1] - length(R);
        R = intersect(S[t], pool);
        m += u;
        push!(num, n[t] * m);
        push!(denom, length(R));
    end
    return sum(num) / sum(denom);
end

function chao(N_o, f)
    return N_o + get(f, 1, 0)^2.0 / (2.0 * get(f, 2, 0));
end

function chao_corrected(N_o, T, f)
    m1 = 2.0 * get(f, 2, 0)  / get(f, 1, 0);
    m2 = 6.0 * get(f, 3, 0) / get(f, 1, 0);
    return (N_o
            + get(f, 1, 0)^2.0 / (2.0 * get(f, 2, 0))
            * (1.0 - m1 / T) / (1.0 - m2 / (T * m1)));
end

function jackknife(N_o, T, f, k)
    if k == 1
        return N_o + (T - 1.0) / T * get(f, 1, 0);
    elseif k == 2
        return (N_o
                + (2.0 * T - 3.0) / T * get(f, 1, 0)
                - (T - 2.0)^2.0 / (T^2.0 - T) * get(f, 2, 0));
    elseif k == 3
        return (N_o
                + (3.0 * T - 6.0) / T * get(f, 1, 0)
                - (3.0 * T^2.0 - 15.0 * T + 19.0) / (T^2.0 - T) * get(f, 2, 0)
                + (T - 3.0)^3.0 / (T * (T - 1.0) * (T - 2.0)) * get(f, 3, 0));
    elseif k == 4
        return (N_o
                + (4.0 * T - 10.0) / T * get(f, 1, 0)
                - (6.0 * T^2.0 - 36.0 * T + 55.0) / (T^2.0 - T) * get(f, 2, 0)
                + (4.0 * T^3.0 - 42.0 * T^2.0 + 148.0 * T - 175.0) / (T * (T - 1.0) * (T - 2.0)) * get(f, 3, 0)
                - (T - 4.0)^4.0 / (T * (T - 1.0) * (T - 2.0) * (T - 3.0)) * get(f, 4, 0));
    elseif k == 5
        return (N_o
                + (5.0 * T - 15.0) / T * get(f, 1, 0)
                - (10.0 * T^2.0 - 70.0 * T + 125.0) / (T^2.0 - T) * get(f, 2, 0)
                + (10.0 * T^3.0 - 120.0 * T^2.0 + 485.0 * T - 660.0) / (T * (T - 1.0) * (T - 2.0)) * get(f, 3, 0)
                - ((T - 4.0)^5.0 - (T - 5.0)^5.0) / (T * (T - 1.0) * (T - 2.0) * (T - 3.0)) * get(f, 4, 0)
                + (T - 5.0)^5.0 / (T * (T - 1.0) * (T - 2.0) * (T - 3.0) * (T - 4.0)) * get(f, 5, 0));
    else
        error("k must be between 1 and 5");
    end
end

function loglh(alpha, N_u, S, O, n, draws)
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
    truncation = 1.0 - monte_carlo(P, comp_P, zeros(T));
    if truncation < 0
        truncation = 5e-200;
    end
    lh = -N_o * log(truncation) + sum_term;
    println("alpha = $alpha, N_u = $N_u, lh = $lh");
    return lh;
end

function simulated_annealing(accepts, sigma_alpha, sigma_Nu,
                             alpha_0, N_u_0, S, O, n, mc_draws)
    alphas = Vector{Float64}();
    N_us = Vector{Float64}();
    evals = Vector{Float64}();
    temps = Vector{Float64}();
    run = 1.0;
    alpha_old = alpha_0;
    N_u_old = N_u_0;
    while length(evals) < accepts
        T = 1.0 / log(1.0 + run);
        lh_old = loglh(alpha_old, N_u_old, S, O, n, mc_draws);
        alpha_new = alpha_old + rand(Normal(0.0, sigma_alpha));
        N_u_new = N_u_old + rand(Normal(0.0, sigma_Nu));
        if (N_u_new < 0)
            N_u_new = 0.0
        end
        if (alpha_new <= 0.0)
            alpha_new = 1e-3
        end
        lh_new = loglh(alpha_new, N_u_new, S, O, n, mc_draws);
        del_lh = lh_new - lh_old;
        rho = minimum([exp(del_lh / T), 1.0]);
        rho = maximum([rho, 0.0]);
        accept = rand(Bernoulli(rho));
        if accept
            alpha_old = alpha_new;
            N_u_old = N_u_new;
            lh_new = loglh(alpha_new, N_u_new, S, O, n, mc_draws);
            push!(alphas, alpha_old);
            push!(N_us, N_u_old);
            push!(evals, lh_new);
            push!(temps, T);
            println("alpha = $alpha_old, N_u = $N_u_old, lh = $lh_new");
        end
        run += 1;
    end
    df = DataFrame(alpha = alphas, N_u = N_us, lh = evals, temp = temps);
    return df;
end

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
