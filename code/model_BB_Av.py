from model_BB import *
from model_BB import ThetaParam as ThetaParamBase


class ThetaParam(ThetaParamBase):
    '''
    Inherit the base model and extend it to include an extra parameter: scattering

    :param grid: parameters corresponding to the dimensions of the grid.
    :type grid: 1D np.array
    '''
    def __init__(self, grid, vz=0.0, vsini=0.0, logOmega=0.0, Av=0.0, T_BB=1100.0, logOmega2=0.0):
        ThetaParamBase.__init__(self, grid, vz, vsini, logOmega, Av, T_BB, logOmega2)

    def __repr__(self):
        return "grid:{} vz:{} vsini:{} logOmega:{} Av:{} T_BB:{} logOmega2:{}".format(self.grid,
         self.vz, self.vsini, self.logOmega, self.Av, self.T_BB, self.logOmega2)