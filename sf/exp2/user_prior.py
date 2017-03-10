def user_defined_lnprior(p):
    '''
    Takes a vector of stellar parameters and returns the ln prior.
    '''
    if not ( (p[2] < 0.5) and (p[2] > -0.5) and
             (p[3] < 1000.0) and (p[3] > -1000.0) and
             (p[4] < 500.0) and (p[4] > 2.0) and
             (p[6] < 1500.0) and (p[6] > 1000.0) and
             (p[11] > 0.0) and (p[11] < 2.0) and
             (p[12] < -1.0) and (p[12] > -2.0) and
             (p[13] > 0.0) and (p[13] < 300 ) ):

    	return -np.inf
    
    # Solar metalicity to within +/- 0.05 dex
    lnp_FeH = -1.0*p[2]**2/(2.0*0.05**2)
    lnp_vz = -1.0*(p[3] - 110.0)**2/(2.0*20.0**2)
    lnp_vsini = -1.0*(p[4] - 105.1)**2/(2.0*30.0**2)
    lnp_l = -1.0*(p[13] - 124.4)**2/(2.0*10.0**2)
    lnp_logAmp = -1.0*(p[12] + 1.40)**2/(2.0*0.03**2)
    lnp_T_BB = -1.0*(p[6] - 1100.0)**2/(2.0*10.0**2)
    # No prior on v_z.
    ln_prior_out = lnp_FeH + lnp_vz + lnp_vsini + lnp_l + lnp_logAmp + lnp_T_BB

    return ln_prior_out

