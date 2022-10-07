using .Utils

using Graphs
using StatsBase
using NLopt

import Random: seed!

NODES = 1000;
EDGES_PER_NODE = 4;
T = [2, 4, 10];
SNOWBALL_WAVES = 2;

function sample_triangles(g, waves, root)
    gs = egonet(g, root, waves);
    s_tri = Set([Set(i) for i in simplecycles_limited_length(gs, 3) if length(i) == 3]);
    return s_tri;
end

function write_row(filename, row)
    open(filename, "a") do io
        for i in row
            if i != row[length(row)]
                print(io, i, ",");
            else
                print(io, i, "\n");
            end
        end
    end
end


G = barabasi_albert(NODES, EDGES_PER_NODE);
estis = [];
for t in T
    for trial in 1:50
        roots = rand(1:NODES, t);
        S = [sample_triangles(G, SNOWBALL_WAVES, r) for r in roots];
        n = [length(s) for s in S];
        O = Set([i for j in S for i in j]);
        LL(x, grad) = -loglh_truncated(x[1], x[2], S, O, t, n, 1000);
        opt = Opt(:LN_SBPLX, 2);
        lower = [0.01, 0];
        opt.lower_bounds = lower;
        opt.min_objective = LL;
        opt.xtol_abs = 0.1;
        (minf, minx, ret) = NLopt.optimize(opt, [5.0, 10.0]);
        push!(estis, [minx[1] ,round(minx[2]) + length(O), t]);
        write_row("tris.csv", [minx[1] ,round(minx[2]) + length(O), t]);
    end
end
