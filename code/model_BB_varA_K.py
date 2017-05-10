from model_BB import *
from model_BB import ThetaParam as ThetaParamBase


class ThetaParam(ThetaParamBase):
    '''
    Inherit the base model and extend it to include an extra parameter: scattering

    :param grid: parameters corresponding to the dimensions of the grid.
    :type grid: 1D np.array
    '''
    def __init__(self, grid, vz=0.0, vsini=0.0, logOmega=0.0, A_K=0.0, T_BB=1100.0, logOmega2=0.0, exponentA_K=-1.999):
        self.grid = grid
        self.vz = vz
        self.vsini = vsini
        self.logOmega = logOmega #log10Omega
        self.T_BB = T_BB
        self.logOmega2 = logOmega2 #log10Omega2
        #Ignore Av-- we are using A_K instead
        self.exponentA_K = exponentA_K
        self.A_K = A_K

    def __repr__(self):
        return "grid:{} vz:{} vsini:{} logOmega:{} A_K:{} T_BB:{} logOmega2:{} exponentA_K:{}".format(self.grid,
         self.vz, self.vsini, self.logOmega, self.A_K, self.T_BB, self.logOmega2, self.exponentA_K)

class ThetaParamDelR(ThetaParamBase):
    '''
    Inherit the base model and extend it to include an extra parameter: scattering

    :param grid: parameters corresponding to the dimensions of the grid.
    :type grid: 1D np.array
    '''
    def __init__(self, grid, vz=0.0, delR=0.0, logOmega=0.0, A_K=0.0, T_BB=1100.0, logOmega2=0.0, exponentA_K=-1.999):
        self.grid = grid
        self.vz = vz
        self.delR = delR
        self.logOmega = logOmega #log10Omega
        self.T_BB = T_BB
        self.logOmega2 = logOmega2 #log10Omega2
        #Ignore Av-- we are using A_K instead
        self.exponentA_K = exponentA_K
        self.A_K = A_K

    def __repr__(self):
        return "grid:{} vz:{} delR:{} logOmega:{} A_K:{} T_BB:{} logOmega2:{} exponentA_K:{}".format(self.grid,
         self.vz, self.delR, self.logOmega, self.A_K, self.T_BB, self.logOmega2, self.exponentA_K)        