#!/usr/bin/env python
import time
import numpy as np
import argparse
parser = argparse.ArgumentParser(prog="star_protostars_varA_K.py", description="Run Starfish fitting model in single order mode with many walkers, with extra black body component, and variable A_K exponent")
parser.add_argument("--samples", type=int, default=5, help="How many samples to run?")
parser.add_argument("--incremental_save", type=int, default=100, help="How often to save incremental progress of MCMC samples.")
parser.add_argument("--resume", action="store_true", help="Continue from the last sample. If this is left off, the chain will start from your initial guess specified in config.yaml.")
args = parser.parse_args()

import os

import multiprocessing
import Starfish.grid_tools
from Starfish.spectrum import DataSpectrum, Mask, ChebyshevSpectrum
from Starfish.emulator import Emulator
import Starfish.constants as C
from Starfish.covariance import get_dense_C, make_k_func, make_k_func_region
from numpy.polynomial import Chebyshev as Ch

from scipy.special import j1
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.linalg import cho_factor, cho_solve
from numpy.linalg import slogdet
from astropy.stats import sigma_clip
from astropy.analytic_functions import blackbody_lambda
from astropy import units as u

import gc
import logging

from itertools import chain
#from collections import deque
from operator import itemgetter
import yaml
import shutil
import json

from star_base import Order as OrderBase
from star_base import SampleThetaPhi as SampleThetaPhiBase 
from model_BB_varA_K import ThetaParam, PhiParam

Starfish.routdir = ""

# list of keys from 0 to (norders - 1)
order_keys = np.arange(1)
DataSpectra = [DataSpectrum.open(os.path.expandvars(file), orders=Starfish.data["orders"]) for file in Starfish.data["files"]]
# list of keys from 0 to (nspectra - 1) Used for indexing purposes.
spectra_keys = np.arange(len(DataSpectra))

#Instruments are provided as one per dataset
Instruments = [eval("Starfish.grid_tools." + inst)() for inst in Starfish.data["instruments"]]

logging.basicConfig(format="%(asctime)s - %(levelname)s - %(name)s -  %(message)s", filename="{}log.log".format(
    Starfish.routdir), level=logging.DEBUG, filemode="w", datefmt='%m/%d/%Y %I:%M:%S %p')

class Order(OrderBase):

    def initialize(self, key):
        OrderBase.initialize(self, key)
        self.extinction = None
        #self.sca_lambda = None


    def update_Theta(self, p):
        OrderBase.update_Theta(self, p)
        # Add in the veiling and scattering
        A_lambda_A_K =  (self.wl/22000.0)**p.exponentA_K #std law is -1.74, Nishiyama is -2.0
        total_extinction = 10**(-0.4 * p.A_K * A_lambda_A_K)
        self.extinction = total_extinction / np.mean(total_extinction) 
        #self.sca_lambda = p.sca*(self.wl/10000.0/2.2)**-4.0 #Rayleigh Scattering


    def evaluate(self):
        '''
        Return the lnprob using the current version of the C_GP matrix, data matrix,
        and other intermediate products.
        '''

        self.lnprob_last = self.lnprob

        X = (self.extinction * self.chebyshevSpectrum.k * self.flux_std * np.eye(self.ndata)).dot(self.eigenspectra.T)

        CC = self.Omega**2 * self.flux_scalar**2 * X.dot(self.C_GP.dot(X.T)) + self.data_mat

        try:
            factor, flag = cho_factor(CC)
        except np.linalg.linalg.LinAlgError:
            print("Spectrum:", self.spectrum_id, "Order:", self.order)
            self.CC_debugger(CC)
            raise

        try:
            # Stellar photosphere 
            model1 = self.Omega * self.flux_scalar * (self.extinction * self.chebyshevSpectrum.k * self.flux_mean + X.dot(self.mus))
            # Black body
            model2 = self.Omega2 * self.BB_lam * self.extinction * self.chebyshevSpectrum.k

            # Return a "metadata blob" of both spectrum models
            raw_models = np.array([model1, model2, self.extinction, self.chebyshevSpectrum.k])

            net_model = model1 + model2
            R = self.fl - net_model

            logdet = np.sum(2 * np.log((np.diag(factor))))
            self.lnprob = -0.5 * (np.dot(R, cho_solve((factor, flag), R)) + logdet)

            self.logger.debug("Evaluating lnprob={}".format(self.lnprob))
            return self.lnprob, raw_models

        # To give us some debugging information about what went wrong.
        except np.linalg.linalg.LinAlgError:
            print("Spectrum:", self.spectrum_id, "Order:", self.order)
            raise


class SampleThetaPhi(Order, SampleThetaPhiBase):
    pass


# Run the program.

model = SampleThetaPhi(debug=True)

model.initialize((0,0))

def lnlike(p):
    # Now we can proceed with the model
    try:
        pars1 = ThetaParam(grid=p[0:3], vz=p[3], vsini=p[4], logOmega=p[5], T_BB=p[6],
                            logOmega2=p[7], A_K=p[8], exponentA_K=p[9])
        model.update_Theta(pars1)
        # hard code npoly=3 (for fixc0 = True with npoly=4)
        pars2 = PhiParam(0, 0, True, p[10:13], p[13], p[14], p[15])
        model.update_Phi(pars2)
        lnp, raw_models = model.evaluate()
        return lnp, raw_models
    except C.ModelError:
        model.logger.debug("ModelError in stellar parameters, sending back -np.inf {}".format(p))
        return -np.inf, []


def lnprior(p):
    print(p[13])
    if not ( p[13] > 0):
        return -1.0*np.inf


# Try to load a user-defined prior
try:
    sourcepath_env = Starfish.config['Theta_priors']
    sourcepath = os.path.expandvars(sourcepath_env)
    with open(sourcepath, 'r') as f:
        sourcecode = f.read()
    code = compile(sourcecode, sourcepath, 'exec')
    exec(code)
    lnprior = user_defined_lnprior
    print("Using the user defined prior in {}".format(sourcepath_env))
except:
    raise


x_vec = np.arange(-1, 1, 0.01)
def cheb_prior(p):
    ch_tot = Ch([0, p[10], p[11], p[12]])
    ch_spec = ch_tot(x_vec)
    if not ( (np.max(ch_spec) < 0.01) and
             (np.min(ch_spec) > -0.01) ):
        return -np.inf

    return 0.0



# Insert the prior here
def lnprob(p):
    #print(p)
    lp0 = lnprior(p)
    lp = lp0 + cheb_prior(p)
    if not np.isfinite(lp):
        return -1.0*np.inf, []
    lnlike_val, raw_models = lnlike(p)
    return lp + lnlike_val, raw_models


import emcee

start = Starfish.config["Theta"]
fname = Starfish.specfmt.format(model.spectrum_id, model.order) + "phi.json"
phi0 = PhiParam.load(fname)

ndim, nwalkers = 16, 40

p0 = np.array(start["grid"] + [start["vz"], start["vsini"], start["logOmega"], start["T_BB"], 
              start["logOmega2"], start["A_K"] , start["exponentA_K"]] +
             phi0.cheb.tolist() + [phi0.sigAmp, phi0.logAmp, phi0.l])

p0_std = [5, 0.02, 0.0005, 0.5, 0.5, -0.01, 5, -0.01, 0.1, 0.0001, -0.005, -0.005, -0.005, 0.01, 0.001, 0.5]

if args.resume:
    p0_ball = np.load("emcee_chain.npy")[:,-1,:]
else:
    p0_ball = emcee.utils.sample_ball(p0, p0_std, size=nwalkers)

n_threads = multiprocessing.cpu_count()
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, threads=n_threads)


nsteps = args.samples
ninc = args.incremental_save
for i, (pos, lnp, state, raw_models) in enumerate(sampler.sample(p0_ball, iterations=nsteps)):
    if (i+1) % ninc == 0:
        time.ctime() 
        t_out = time.strftime('%Y %b %d,%l:%M %p') 
        print("{0}: {1:}/{2:} = {3:.1f}%".format(t_out, i, nsteps, 100 * float(i) / nsteps))
        
        np.save('temp_emcee_chain.npy',sampler.chain)
        
        try:
            np.save('temp_raw_models.npy',sampler.blobs)
        except:
            print("The temporary raw models could not save, something is off.")

np.save('emcee_chain.npy',sampler.chain)
np.save('raw_models.npy',sampler.blobs)

print("The end.")