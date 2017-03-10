def user_defined_lnprior(p):
    '''
    Takes a vector of stellar parameters and returns the ln prior.
    '''
    if not ( (p[12] > 0.0) and 
             (p[10] > 0.0) and
             (p[3] < 1000.0) and (p[3] > -1000.0) and
             (p[4] < 500.0) and (p[4] > 2.0)):
    	return -np.inf
    
    # Solar metalicity to within +/- 0.05 dex
    lnp_FeH = -1.0*p[2]**2/(2.0*0.05**2)
    # v_z=0 +/- 100 km/s, made up to improve convergence
    lnp_vz = -1.0*(p[3] - 0.0)**2/(2.0*100.0**2)
    # No prior on v_z.
    ln_prior_out = lnp_FeH + lnp_vz

    return ln_prior_out

