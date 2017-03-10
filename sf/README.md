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
