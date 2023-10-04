struct HyperParams
    F2::uwreal
    F3::uwreal 
    Lambda::uwreal
end
function HyperParams(; _F2=uwreal([1.7505,0.0089], "F2"), F3_val=0.5226, Lambda_val=341, Lambda_F3_cov= [[150.7862014579557, -0.01716315895017156] [ -0.01716315895017156, 0.00001830314004469306]])
    vals = cobs([Lambda_val, F3_val], Lambda_F3_cov, [1,2])
    return HyperParams(_F2, vals[2], vals[1])
end

"""
beta_function_coeff(nc, nf; nl=5)
Coefficient of the beta function at nl=4,5 loops (default is nl=5) with nc colors and nf flavours 
"""
function beta_function_coeff(nc, nf; nl=5)
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

    if nl ==5
        push!(bcoefs, ( 524.56 - 181.8 * nf + 17.16 * nf^2 - 0.22586 * nf^3 - 0.0017993 * nf^4 ) / (4 * pi^2)^5 )
    end
    return bcoefs
end


"""
tau_function_coeff(nc, nf; nl=5)
Coefficient of the tau function at nl=4,5 loops (default is nl=5) with nc colors and nf flavours 
"""
function tau_function_coef(nc, nf; nl=5)
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

    if nl == 5
        push!(tcoefs, 2*4^(4+1)/(4*pi)^(2*(4+1)) * (559.7069 - 143.6864*nf + 7.4824*nf^2 + 0.1083*nf^3 - 0.000085359*nf^4) )
    end
    
    return tcoefs
end
