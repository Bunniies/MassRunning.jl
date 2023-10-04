

"""
error_mass_running(obs::uwreal; F4=uwreal([0.6824, 0.0052], "F4"); hp::HyperParams=HyperParams())

Given an observable in the SF scheme at scale muhad and the running factor F4, this function performs
linear error propagation
"""
function error_mass_running(obs::uwreal; F4=uwreal([0.6824, 0.0052], "F4"), hp::HyperParams=HyperParams())

    F2 = hp.F2
    F3 = hp.F3 
    Lamb = hp.Lambda  
    Cov_elem_12 =  -0.01716315895017156
    uwerr(obs); uwerr(F2); uwerr(F3); uwerr(F4); uwerr(Lamb)

    dF4dLambda = err(F4) / err(Lamb)
    
    delta_obs2 = (
                err(obs)^2 * value(F2)^2 * value(F3)^2 * value(F4)^2 +
                value(obs)^2 * err(F2)^2 * value(F3)^2 * value(F4)^2 +
                value(obs)^2 * value(F2)^2 * err(F3)^2 * value(F4)^2 +
                value(obs)^2 * value(F2)^2 * value(F3)^2 * err(F4)^2 +
                2 * value(obs)^2 * value(F2)^2 * value(F3) * value(F4) *  Cov_elem_12  * dF4dLambda
    )    
    return delta_obs2^0.5
end

"""
msbar_over_MRGI_factor(; mu=uwreal([3000,0.0],"mu"), nc = 3, nf = 4, hp::HyperParams=HyperParams())

Given a scale mu. Given nc number of color and nf number of flavours for the beta and tau functions. Given Lambda included in hp.

This function returns the F4 coefficient required to perform the RGI running
"""
function msbar_over_MRGI_factor(; mu=uwreal([3000,0.0],"mu"), nc=3, nf=4, nl=5, hp::HyperParams=HyperParams())

    bcoefs = beta_function_coeff(nc, nf, nl=nl)
    tcoefs = tau_function_coef(nc, nf, nl=nl)
    
    gbar = g_from_RG_eq(mu, bcoefs, hp=hp, nl=nl)
    F4 = solve_RG_eq_mass(gbar, bcoefs, tcoefs, nl=nl)
    # uwerr(gbar); uwerr(F4)
    # println("g  = ", gbar)
    # println("F4 = ", F4)
    return F4
end


function solve_RG_eq(g, bt; nl=5)
    if nl == 5
        beta = (x,p) -> - x^3*p[1] - x^5*p[2] - x^7*p[3] - x^9*p[4] - x^11*p[5] 
    else
        beta = (x,p) -> - x^3*p[1] - x^5*p[2] - x^7*p[3] - x^9*p[4] 
    end

    integrand(x,p) =  1 / beta(x,p) + 1 / ( p[1]*x^3) - p[2]/(p[1]^2*x)
    N = (bt[1] * g^2)^(bt[2]/(2*bt[1]^2)) * exp(1 / (2*bt[1]*g^2))
    int_val = int_error(integrand, 0.0, g, uwreal.(bt) )
    return N * exp(int_val)
end

function g_from_RG_eq(mu, bcoef; hp::HyperParams=HyperParams(), nl=5, grange=1.5:0.001:2.5, pl::Bool=false)
    mu_over_lambda_predict = []
    grange1 =[]
    mu_over_lambda_true = mu / hp.Lambda ; uwerr(mu_over_lambda_true)

    for g in grange
        try
            aux = solve_RG_eq(g, bcoef, nl=nl)
            push!(mu_over_lambda_predict, aux)
            push!(grange1, g)
        catch
        end
    end
    
    @. model(x, p)  = p[1] + p[2]*x  + p[3]*x^2 + p[4]*x^3 + p[5]* x^5 + p[6] * x^6 + p[7] * x^7 + p[8] * x^8 
    fit = curve_fit(model, grange1, value.(mu_over_lambda_predict), fill(0.5,8))
    f(x,p)=  p[1] .- model(x,coef(fit))
    gbar = root_error(f, 2.0, [mu_over_lambda_true] ) ; uwerr(gbar)
    if pl
        plot(grange1, model(grange1, coef(fit)))
        v = value(mu_over_lambda_true)
        e = err(mu_over_lambda_true)
        plot(grange1, value.(mu_over_lambda_predict))
        errorbar(value(gbar), xerr=err(gbar),  value(mu_over_lambda_true), yerr=err(mu_over_lambda_true), fmt="d", color="black", capsize=2)
        fill_between(grange1, v-e, v+e, alpha=0.5, color="tomato")
        ylabel(L"$\mu/\Lambda$")
        xlim(value(gbar) - 0.2, value(gbar) + 0.2)
        ylim(value(mu_over_lambda_true) - 3, value(mu_over_lambda_true) +3)
        xlabel(L"$g$")
        display(gcf())
        close("all")
    end

    return gbar
end

function solve_RG_eq_mass(g, bt, dt; nl=5)
    if nl == 5 
        beta = (x,p) -> -x^3*p[1] - x^5*p[2] - x^7*p[3] - x^9*p[4] - x^11*p[5]
        tau = (x,p) -> -x^2*p[6] - x^4*p[7] - x^6*p[8] - x^8*p[9] - x^10*p[10]
    else
        beta = (x,p) -> -x^3*p[1] - x^5*p[2] - x^7*p[3] - x^9*p[4] 
        tau = (x,p) -> -x^2*p[5] - x^4*p[6] - x^6*p[7] - x^8*p[8] 
    end

    integrand(x,p) = tau(x,p) / beta(x,p) - (dt[1]/(bt[1]*x)) 
    int_val = int_error(integrand, 1e-10, g, uwreal.(vcat(bt, dt)) )
    N = (2*bt[1]*g^2)^(dt[1]/(2*bt[1]))

    return N * exp(int_val)
end

"""
running_to_scale_invariant_mass(guess::uwreal, Mrgi::uwreal; nc=3, nf=3, nl=5, hp::HyperParams=HyperParams())

Starting from an intial guess and a RGI mass, this functions performs the running to the scale invariant mass m_q(mu=m_q, N_f=3).
The beta and tau function coeffcients are computed with nc and nf numbers of colors and flavours, respectively.
The perturbative expressions for the beta and tau functions are solved with nl number of loops.
"""
function running_to_scale_invariant_mass(guess::uwreal, Mrgi::uwreal; nc=3, nf=3, nl=5, hp::HyperParams=HyperParams())
    F4 = msbar_over_MRGI_factor(mu=guess, nc=nc, nf=nf, hp=hp)
    msi_prec = Mrgi * F4
    while true
        F4 = msbar_over_MRGI_factor(mu=msi_prec, nc=nc, nf=nf, hp=hp)
        msi = Mrgi * F4
        if abs(value((msi-msi_prec)/msi)) < 0.0001
            return msi
        end 
        msi_prec = msi
    end
end
"""
fzeta(r, r0)

Matching between N_f=3 -> N_f=4
Parametrization of the decoupling function of the mass: zeta[mu=msi3](r) = m[MSbar,mu=mc_si3,N_f=4] / mbar[MSbar,mu=m_si3,N_f=3]
    ## zeta(r) = zeta03 + zeta13 * (r-r0) + zeta23 * (r-r0)^2
"""
function fzeta(r, r0)
    zeta03 =  0.992628
    zeta13 =  0.00344098
    zeta23 = -0.00165175
    tmp = zeta03 + zeta13*(r-r0) + zeta23*(r-r0)^2
    return tmp
end

