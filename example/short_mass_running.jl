using ADerrors

function tau_over_beta(g2_1, g2_2)
    CC_f = [[2.33798e-2, -1.47011e-2, 2.81966e-3, -1.66404e-4] [ -1.47011e-2, 9.54563e-3, -1.87752e-3, 1.12962e-4] [2.81966e-3, -1.87752e-3, 3.78680e-4, -2.32927e-5] [-1.66404e-4, 1.12962e-4, -2.32927e-5, 1.46553e-6] ]
    p = cobs([1.28493, -0.292465, 0.0606401, -0.00291921], CC_f, "tau_over_beta")
    func(x,p) = 1 / x * (p[1] + p[2]*x^2+ p[3]x^4 + p[4]x^6)
    res = int_error(func, g2_1, g2_2, p)
    return exp(res)

end

g1 = uwreal([3.949, 0.011 ], "g0")
# g1 = uwreal([3.949, 0.000 ], "g0")

g2 = uwreal([2.6723, 0.0064], "g_mu0/2")
F3 = tau_over_beta(sqrt(g1), sqrt(g2)); uwerr(F3); F3
F2 = uwreal([1.7505, 0.0089], "M/m_mu0")
rgi_factor = F2 * F3 ; uwerr(rgi_factor); rgi_factor