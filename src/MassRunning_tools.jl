struct HyperParams
    F2::uwreal
    F3::uwreal 
    Lambda::uwreal
    match_nf4_to_nf3::Float64 
end
function HyperParams(; _F2=uwreal([1.7505,0.0089], "F2"), F3_val=0.5226, Lambda_val=341, Lambda_F3_cov= [[150.7862014579557, -0.01716315895017156] [ -0.01716315895017156, 0.00001830314004469306]], nf4_to_nf3=0.9929)
    vals = cobs([Lambda_val, F3_val], Lambda_F3_cov, [1,2])
    return HyperParams(_F2, vals[2], vals[1], nf4_to_nf3)
end

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
function msbar_over_MRGI_factor(; mu=uwreal([3000,0.0],"mu"), nc = 3, nf = 4, hp::HyperParams=HyperParams())

    bcoefs = beta_function_coeff(nc, nf)
    tcoefs = tau_function_coef(nc, nf)
    
    gbar = g_from_RG_eq(mu, bcoefs, hp=hp)
    F4 = solve_RG_eq_mass(gbar, bcoefs, tcoefs)
    # uwerr(gbar); uwerr(F4)
    # println("g  = ", gbar)
    # println("F4 = ", F4)
    return F4
end


function solve_RG_eq(g, bt)
    beta(x,p) = - x^3*p[1] - x^5*p[2] - x^7*p[3] - x^9*p[4] 
    integrand(x,p) =  1 / beta(x,p) + 1 / ( p[1]*x^3) - p[2]/(p[1]^2*x)
    N = (bt[1] * g^2)^(bt[2]/(2*bt[1]^2)) * exp(1 / (2*bt[1]*g^2))
    int_val = int_error(integrand, 0.0, g, uwreal.(bt) )
    return N * exp(int_val)
end

function g_from_RG_eq(mu, bcoef;hp::HyperParams=HyperParams(), grange=1.5:0.001:2.5, pl::Bool=true)
    mu_over_lambda_predict = []
    grange1 =[]
    mu_over_lambda_true = mu / hp.Lambda ; uwerr(mu_over_lambda_true)

    for g in grange
        try
            aux = solve_RG_eq(g, bcoef)
            push!(mu_over_lambda_predict, aux)
            push!(grange1, g)
        catch
        end
    end
    
    @. model(x, p)  = p[1] + p[2]*x  + p[3]*x^2 + p[4]*x^3 + p[5]* x^5 + p[6] * x^6 + p[7] * x^7
    fit = curve_fit(model, grange1, value.(mu_over_lambda_predict), fill(0.5,7))
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

function solve_RG_eq_mass(g, bt, dt)
    beta(x,p) = -x^3*p[1] - x^5*p[2] - x^7*p[3] - x^9*p[4] #- x^11*p[5]
    tau(x, p) = -x^2*p[5] - x^4*p[6] - x^6*p[7] - x^8*p[8] #- x^10*p[10]
    integrand(x,p) = tau(x,p) / beta(x,p) - (dt[1]/(bt[1]*x)) 

    int_val = int_error(integrand, 1e-20, g, uwreal.(vcat(bt, dt)) )
    N = (2*bt[1]*g^2)^(dt[1]/(2*bt[1]))

    return N * exp(int_val)

end

"""
beta_function_coeff(nc, nf)
Coefficient of the beta function at 4 loops with nc colors and nf flavours 
"""
function beta_function_coeff(nc, nf)
    zeta3 = 1.20205690315959428539; zeta4 = pi^4/90;  zeta5 = 1.03692775514336992633;

    tf = 1/2; ca = nc; cf = (nc^2 - 1)/(2*nc); na = nc^2 - 1; nr = nc; 
    dada = nc^2/24*(nc^2 + 36)*na;
    dfda = (1/48)*(nc - 1)*(nc + 1)*(nc^2 + 6)*nr;
    dfdf = (nc - 1)*(nc + 1)*(nc^4 - 6*nc^2 + 18)/96/nc^3*nr; 

    bcoefs = [((11/3*ca - 4/3*tf*nf)/(4*pi)^2), 
    (1/(4*pi)^4)*(34/3*ca^2 - 4*cf*tf*nf - 20/3*ca*tf*nf), 
    (1/(4*pi)^6)*(2857/54*ca^3 + 2*cf^2*tf*nf - 205/9*cf*ca*tf*nf - 1415/27*ca^2*tf*nf + 44/9*cf*tf^2*nf^2 + 158/27*ca*tf^2*nf^2), 
    (1/(4*pi)^8)*(ca^4*(150653/486 - 44/9*zeta3) + ca^3*tf*nf*(-39143/81 + 136/3*zeta3) + ca^2*cf*tf*nf*(7073/243 - 656/9*zeta3) +ca*cf^2*tf*nf*(-4204/27 + 352/9*zeta3) + 46*cf^3*tf*nf +
    ca^2*tf^2*nf^2*(7930/81 + 224/9*zeta3) + cf^2*tf^2*nf^2*(1352/27 - 704/9*zeta3) + ca*cf*tf^2*nf^2*(17152/243 + 448/9*zeta3) + 424/243*ca*tf^3*nf^3 +
    1232/243*cf*tf^3*nf^3 + (dada/na)*(-80/9 + 704/3*zeta3) + nf*(dfda/na)*(512/9 - 1664/3*zeta3) +  nf^2*(dfdf/na)*(-704/9 + 512/3*zeta3))] 

    return bcoefs
end

"""
tau_function_coeff(nc, nf)
Coefficient of the tau function at 4 loops with nc colors and nf flavours 
"""
function tau_function_coef(nc, nf)
    zeta3 = 1.20205690315959428539; zeta4 = pi^4/90;  zeta5 = 1.03692775514336992633;

    tf = 1/2; ca = nc; cf = (nc^2 - 1)/(2*nc); na = nc^2 - 1; nr = nc; 
    dfda = (1/48)*(nc - 1)*(nc + 1)*(nc^2 + 6)*nr;
    dfdf = (nc - 1)*(nc + 1)*(nc^4 - 6*nc^2 + 18)/96/nc^3*nr; 

    tcoefs = [(1/(4*pi)^2)*6*cf, 
    (1/(4*pi)^4)*(3*cf^2 + 97*cf*ca/3 - 20*cf*tf*nf/3),
    (1/(4*pi)^6)*(129*cf^3 - 129*cf^2*ca/2 + 11413*cf*ca^2/54 + cf^2*tf*nf*(-92 + 96*zeta3) + cf*ca*tf*nf*(-1112/27 - 96*zeta3) - 280*cf*tf^2*nf^2/27),
    (1/(4*pi)^8)*(cf^4*(-1261/4 - 672*zeta3) +cf^3*ca*(15349/6 + 632*zeta3) +cf^2*ca^2*(-34045/18 - 304*zeta3 + 880*zeta5) +
    cf*ca^3*(70055/36 + 2836*zeta3/9 - 880*zeta5) +cf^3*tf*nf*(-560/3 + 1104*zeta3 - 960*zeta5) +cf^2*ca*tf*nf*(-17638/27 + 736*zeta3 - 528*zeta4 + 160*zeta5) +
    cf*ca^2*tf* nf*(-65459/81 - 5368*zeta3/3 + 528*zeta4 + 800*zeta5) +cf^2*tf^2*nf^2*(608/27 - 320*zeta3 + 192*zeta4) +cf*ca*tf^2*nf^2*(2684/81 + 320*zeta3 - 192*zeta4) +
    cf*tf^3* nf^3*(-1328/81 + 256*zeta3/9) + (dfda/nr)*(-64 + 480*zeta3) + nf*(dfdf/nr)*(128 - 960*zeta3))]

    return tcoefs
end