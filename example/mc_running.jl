using Revise, MassRunning, ADerrors

HyperParams()

##############
## Mc
################
# mc at muhad
mCmuhad = uwreal([1627.0,8.0],"mCmuhad")

mCRGI = uwreal([1455.0, 20.0], "RGI")
# matching nf=3 to nf = 4

match_nf = 1.29011 / 1.29931 # nf4 / nf3
##################
## mc at 3GeV nf=3 the matched 
##################
# test with full running 
mc_3gev_3f = match_nf * mbar_running(mCmuhad, mu=uwreal([3000.0,0.0],"mu"), nf=3); uwerr(mc_3gev_3f);# mc_3gev_3f
# test with error propagation
err_mc_3gev_3f = match_nf *  error_mass_running(mCmuhad, F4=F4=uwreal([0.6824, 0.0052], "F4")) #; uwerr(test); details(test)

println("From full running: ", "m_c(3GeV, N_f=4)= "  , mc_3gev_3f )
println("From linear prop: "," δm_c(3GeV, N_f=3)= "  , err_mc_3gev_3f )
##
#####################
## mc at mc nf=3 then matched 
#####################
# test with full running 
mc_mc_3f  =  match_nf *   mbar_running(mCmuhad, mu=uwreal([1299.9,0.0], "mu"), nf=3); uwerr(mc_mc_3f); mc_mc_3f
# test with error propagation
err_mc_mc_3f =  match_nf *error_mass_running(mCmuhad, F4=uwreal([0.874, 0.013], "F4")) #; uwerr(test); details(test)

println("From full running: ", "m_c(m_c, N_f=3)= ", mc_mc_3f )
println("From linear prop: "," δm_c(m_c, N_f=3)= "  ,err_mc_mc_3f )


##############
## Mb
##############
mBmuhad = uwreal([6604.0,59.0],"mCmuhad")
mCRGI = uwreal([1455.0, 20.0], "RGI")
mb_3gev_3f = match_nf * mbar_running(mBmuhad, mu=uwreal([3000.0,0.0],"mu"), nf=3); uwerr(mb_3gev_3f); mb_3gev_3f
mb_mb_3f   = match_nf *  mbar_running(mBmuhad, mu=uwreal([4170.4,0.0], "mu"), nf=3); uwerr(mb_mb_3f); mb_mb_3f

details(mb_mb_3f)

sqrt(err(mb_mb_3f)^2 * 0.3515)
sqrt(err(mb_mb_3f)^2 * 0.6485)

sqrt(err(mb_3gev_3f)^2 - err(mb_3gev_3f)^2 * 0.57) 