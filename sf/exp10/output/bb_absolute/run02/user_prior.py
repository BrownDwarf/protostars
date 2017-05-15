def user_defined_lnprior(p):
    '''
    Takes a vector of stellar parameters and returns the ln prior.
    '''
    if not ( (p[2] < 0.5) and (p[2] > 0) and
             (p[3] < 1000.0) and (p[3] > -1000.0) and
             (p[4] < 10000.0) and (p[4] > 0.01) and
             (p[6] < 1500.0) and (p[6] > 1000.0) and
             (p[8] < 20.0) and (p[8] > 0.0) and
             (p[9] < -1.7) and (p[9] > -2) and
             (p[13] > 0.0) and (p[13] < 2.0) and
             (p[14] < -1.0) and (p[14] > -3.0) and
             (p[15] > 0.0) and (p[15] < 300 ) ):

    	return -np.inf
    
    # Solar metalicity to within +/- 0.05 dex
    lnp_FeH = -1.0*(p[2] - 0.01)**2/(2.0*0.001**2)
    lnp_vz = -1.0*(p[3] - 110.0)**2/(2.0*20.0**2)
    lnp_l = -1.0*(p[15] - 124.4)**2/(2.0*10.0**2)
    lnp_logAmp = -1.0*(p[14] + 1.40)**2/(2.0*0.03**2)
    # No prior on v_z.
    ln_prior_out = lnp_FeH + lnp_vz + lnp_l + lnp_logAmp

    return ln_prior_out

