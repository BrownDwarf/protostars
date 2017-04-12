## Starfish

This directory contains the configuration files and, in some cases, the MCMC posterior chains for `Starfish`.  This meaning of the files can be found by reading the [Starfish documentation](http://iancze.github.io/Starfish/current/index.html).

```bash
.
├── exp1
│   ├── libraries
│   ├── output
│   │   └── veil1
│   │       └── run01
                ├── config.yaml         <-- Starfish config file
                ├── emcee_chain.npy     <-- 5000 x 40 MCMC chain
└── exp2
    ├── libraries
    ├── output
    │   └── veil1
    │       ├── run01
    │       │   ├── config.yaml
    │       │   └── temp_emcee_chain.npy   <-- failed 1500 x 40 MCMC
    │       └── run02
    │           ├── config.yaml
    │           └── emcee_chain.npy        <-- 5000 x 40 MCMC chain
    │           └── temp_raw_models.npy    <-- all 5000 x 40 models!
    └── plots
```

### First attempt at veiling: *exp1*

We experimented with the `star_veil.py` code, a pre-existing code written with Ben Kidder at McDonald observatory.  This code models the veiling as a constant in `f_lambda` units.  The 5000 step MCMC chain converged, and seems to be a "good" representation of the data.  
However, we abandoned this parameterization in favor of a Black Body with a solid angle model that offers more physical interpretability.

### Second attempt at veiling: *exp2*
Here we parameterize veiling as described in [Starfish Issue #75](https://github.com/iancze/Starfish/issues/75), namely a disk temperature `T_BB` and a disk solid angle, `Omega_BB`.  MGS forked `star_veil.py` into `star_BB.py` to adjust for these new parameters.  
**run01** failed due to a bug in which only the stellar model was being multiplied by the Chebyshev polynomials, leading to highly degenerate posteriors and numerical precision problems that crashed the code 1500 steps into the 5000 step MCMC sampling.  
**run02** corrected the Chebyshev bug and applied some additional strong posteriors to aid in convergence.  These priors are probably overly strong, and should be relaxed in future runs.  The samples converged.

### Third (and fifth) attempts at veiling: *exp3*
Veiling is parameterized as disk temperature `T_BB` and a disk solid angle, `Omega_BB`but with corrected absolute flux ratio using the "flux_scalars" discussion in [Issue 5](https://github.com/BrownDwarf/protostars/issues/5).   
**run01** The whole available spectral region, with a tight T_BB prior of 1100 K.  
**run02** The whole available spectral region, with a relaxed T_BB prior 1000-1700 K.  

### Fourth attempt at veiling: *exp4*
Same as *exp3*, but with a reduced spectral chunk (or "window") surrounding the CO bandhead.  
**run01** A tight T_BB prior of 1100 K, priors the same as exp3 run01  

### Sixth attempt at veiling: *exp5*
Same as *exp3*, but with a reduced spectral chunk (or "window") surrounding the CO bandhead plus Na I line at 2.2 um.   
**run01** A tight T_BB prior of 1100 K, priors the same as exp3 run01.  

### Seventh attempt at veiling: *exp6*
Same as *exp3* run01, but with extinction and scattering implemented
**run01** A tight T_BB prior of 1100 K, priors the same as exp3 run01.  
