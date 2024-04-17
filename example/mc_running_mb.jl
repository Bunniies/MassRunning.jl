using  MassRunning, ADerrors
# running for the b quark mass at the scale L0

# HyperParameters definition 
# for N_f=3 theory
hp  = HyperParams(; _F2=uwreal([1.7505,0.0089], "F2"), F3_val=0.8423, Lambda_val=341, Lambda_F3_cov=[[150.7862014579557, 0.0] [ 0.0, 6.76e-6]])

# for N_f=4 theory
hp4 = HyperParams(; _F2=uwreal([1.,0.0], "tmp"), F3_val=0.5226, Lambda_val=341*0.874415, Lambda_F3_cov=[[150.7862014579557*0.874415^2, 0.0] [ 0.0, 6.76e-6]])

# reference scale 3 GeV
mucref = 3000.0 # MeV

####################################################################
# Inputs
####################################################################

#for charm (start from mSF)
mSFmuhad3 = uwreal([4481.8,20.6],"mSFmuhadNf3") # input: m_q[SF, mu=mu_had, N_f=3]  
M3 = mSFmuhad3 * hp.F2 * hp.F3; uwerr(M3); M3 # RGI mass N_f=3



#####################################################################
# Running from M[N_f=3] to scale invariant ("si") mass, m_q(mu=m_q, N_f=3), through an iterative procedure:
# for  charm
# msi3 = running_to_scale_invariant_mass(uwreal([1297.592027364731,0.0],"mu"), M3, nc=3, nf=3, nl=5, hp=hp); uwerr(msi3) # m_q(mu=m_q[N_f=3], N_f=3)
# for bottom
msi3 = running_to_scale_invariant_mass(uwreal([4198.,0.0],"mu"), M3, nc=3, nf=3, nl=5, hp=hp); uwerr(msi3) # m_q(mu=m_q[N_f=3], N_f=3)
# Matching N_f=3 -> N_f=4
r0 = M3 / hp.Lambda
m4mumc3 = fzeta(r0, value(r0)) * msi3; uwerr(m4mumc3) # m_c(mu=m_c[N_f=3], N_f=4)

println("####################################################")
println("## All masses in MeV")
println("M(RGI, N_f=3) =           ", M3)
println("mc(mu=mc[N_f=3], N_f=3) = ", msi3)
println("mc(mu=mc[N_f=3], N_f=4) = ", m4mumc3)
println("####################################################")


#####################################################################
# Running from M[N_f=3] to m_q(mu=3GeV, N_f=3)
# N_f=3: running RGI --> mu=3 GeV
F4_mu_3 = msbar_over_MRGI_factor(mu=uwreal([mucref,0.0],"mu"), nc=3, nf=3, nl=5, hp=hp)
mmu3 = M3*F4_mu_3; uwerr(mmu3)

# N_f=4: mu=msi3 --> RGI
F4_msi4 = msbar_over_MRGI_factor(mu=msi3, nc=3, nf=4, nl=5, hp=hp4)
M4  = m4mumc3/F4_msi4; uwerr(M4)

# N_f=4: RGI --> mu=3 GeV
F4_mu_4 = msbar_over_MRGI_factor(mu=uwreal([mucref,0.0],"mu"), nc=3, nf=4, nl=5, hp=hp4)
mmu4 = M4*F4_mu_4; uwerr(mmu4)

###############################################################################
# N_f=4: RGI --> mu=m_c[mu=m_c,N_f=4] with iterative procedure
#fo
# for charm
# msi4 = running_to_scale_invariant_mass(uwreal([1297.5,0.0],"mu"), M4, nc=3, nf=4, nl=5, hp=hp4); uwerr(msi4) # m_q(mu=m_q[N_f=4], N_f=4)
# for bottom
msi4 = running_to_scale_invariant_mass(uwreal([4170.,0.0],"mu"), M4, nc=3, nf=4, nl=5, hp=hp4); uwerr(msi4) # m_q(mu=m_q[N_f=4], N_f=4)

#####################################################################

println("####################################################")
println("## mbar[mu] at mu = mu_ref = ", mucref, " MeV\n")
println("mc(mu=mu_ref, N_f=3) =     ", mmu3)
println("M(RGI, N_f=4) =           ", M4)
println("mc(mu=mu_ref GeV, N_f=4) =     ", mmu4)
println("mc(mu=mc[N_f=4], N_f=4) = ", msi4)
println("####################################################")


