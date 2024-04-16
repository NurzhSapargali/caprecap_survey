using Plots

gf(a, N, n) = 1.0 + n / (a * N)
dgf(a, N, n) = -n / (N * a^2)

function f(a)
    g = gf(a, 1, 1)
    dg = dgf(a, 1, 1)
    lg = log(g)
    del_1 = g^a * (lg + a / g * dg) * lg + g^(a - 1) * dg
    del_2 = g^a * (lg + a / g * dg)
    del_3 = g^(a - 1) * (lg + (a - 1) / g * dg)
    del_4 = g^a + a * del_2 - 1
    return ( 1 / a + (del_1 * (g^a - 1) - g^a * lg * del_2) / (g^a - 1)^2
    - (del_3 * (g^a - 1) * a - g^(a - 1) * del_4) / (a^2 * (g^a - 1)^2) )
end