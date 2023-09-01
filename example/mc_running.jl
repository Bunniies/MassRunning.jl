using  MassRunning, ADerrors

# HyperParameters definition 
hp = HyperParams()
hp.F2               # F2 coefficient 
hp.F3               # F3 coefficient 
hp.Lambda           # Lambda from ALPHA (covariance with F3 taken into account)
hp.match_nf4_to_nf3 # matching between nf4 and nf3 theories 

##############
##  Running of Mc
################
# input: mc at muhad
mCmuhad = uwreal([1627.0,8.0],"mCmuhad") # = McRGI * RGI_factor

# running to 3 GeV
F4_3gev = msbar_over_MRGI_factor(mu=uwreal([3000.0,0.0],"mu"), nc=3, nf=3, hp=hp)
mc_3gev_3f = mCmuhad * hp.F2 * hp.F3 * F4_3gev; uwerr(mc_3gev_3f)
mc_3gev_4f = hp.match_nf4_to_nf3 * mc_3gev_3f; uwerr(mc_3gev_4f)
# check with linear error propagation
err_mc_3gev_3f = error_mass_running(mCmuhad, F4=F4_3gev, hp=hp)

println("\n######################## ")
println("Results at 3GeV\n")
println("mc(3GeV, nf=3) = ", mc_3gev_3f)
println("err from lin. prop. =", err_mc_3gev_3f)
println("mc(3GeV, nf=4) = ", mc_3gev_4f)
println("######################## \n")

# running to mass scale mbar
F4_mbar = msbar_over_MRGI_factor(mu=uwreal([1299.0,0.0],"mu"), nc=3, nf=3, hp=hp)
mc_mc_3f = mCmuhad * hp.F2 * hp.F3 * F4_mbar; uwerr(mc_mc_3f)
mc_mc_4f = hp.match_nf4_to_nf3 * mc_mc_3f; uwerr(mc_mc_4f)
# check with linear error propagation
err_mc_mc_3f = error_mass_running(mCmuhad, F4=F4_mbar, hp=hp)

println("\n######################## ")
println("Results at mc scale\n")
println("mc(mc, nf=3) = ", mc_mc_3f)
println("err from lin. prop. =", err_mc_mc_3f)
println("mc(mc, nf=4) = ", mc_mc_4f)
println("########################")

##
"""
You can use the details() function to access the error contribution of each quantity.
Due to a bad behaviour in ADerrors, when introducing covariance between Lambda and F3, the ensemble id is lost in an 
unpredictable integer chosen by aderrors. To overcome this, I had to call the ensembles correspnding to Lamdbda:"1" and the 
ensembles corresponding to F3:"2".
"""
details(mc_3gev_4f)
details(mc_mc_4f)
