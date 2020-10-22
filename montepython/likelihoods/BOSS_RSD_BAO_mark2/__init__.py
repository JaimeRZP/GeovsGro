import os
import numpy as np
from montepython.likelihood_class import Likelihood
import montepython.io_mp as io_mp
import warnings
import math 
import scipy.constants as consts

class BOSS_RSD_BAO_mark2(Likelihood):

    # initialization routine

    def __init__(self, path, data, command_line):

        Likelihood.__init__(self, path, data, command_line)

        #Check for incompatibilities
        if 'kids450_bs_2cosmos_mark2' or 'kv450_bs_2cosmos_mark2' not in data.experiments:
            self.need_cosmo1_arguments(data, {'output': 'mPk'})
            self.need_cosmo1_arguments(data, {'P_k_max_h/Mpc': '1.'})
            self.need_cosmo1_arguments(data, {'z_max_pk': '1.'})

        # define array for values of z and data points
        self.data = np.array([], 'float64')
        self.error = np.array([], 'float64')
        self.z = np.array([], 'float64')
        self.type = np.array([], 'int')
        self.rd_fid = 147.78
        #read covariance matrix
        self.cov_data = np.loadtxt(os.path.join(self.data_directory, self.cov_file))

        # read redshifts and data points
        with open(os.path.join(self.data_directory, self.data_file), 'r') as filein:
            for line in filein:
                if line.strip() and line.find('#') == -1:
                    # the first entry of the line is the identifier
                    this_line = line.split()
                    # insert into array if this id is not manually excluded
                    self.z        = np.append(self.z, float(this_line[1]))
                    self.data     = np.append(self.data, float(this_line[2]))
                    self.type     = np.append(self.type, int(this_line[3]))
                    try:
                        self.error = np.append(self.error, float(this_line[4]))
                        print('error added')
                    except:
                        print('included in covariance')

        return 


    def loglkl(self, cosmo1, cosmo2, data):
        #ZP: Being the position of the acoustic peak 
        # a strong goemetry phenomenon, calculations 
        #are done using cosmo2 as it is the one responsible to 
        #parametrize goemetry. 

        #However, both cosmologies are called for matters 
        #of consistency at the sampler 

        chi2 = 0.

        # for each point, compute angular distance da, radial distance dr,
        # volume distance dv, sound horizon at baryon drag rs_d,
        # theoretical prediction and chi2 contribution
        diff = []
        for counter, item in enumerate(self.data):
            

            if self.type[counter] == 1:
                theo = cosmo1.scale_independent_growth_factor_f(self.z[counter])*cosmo1.sigma(8./cosmo1.h(),self.z[counter]) 
                
            elif self.type[counter] == 2:
                theo = self.Dv(cosmo2, self.z[counter])*self.rd_fid/cosmo2.rs_drag()

            elif self.type[counter] == 3:
                theo = cosmo2.rs_drag()/self.Dv(cosmo2, self.z[counter])

            elif self.type[counter] == 4:
                theo = cosmo2.angular_distance(self.z[counter])/cosmo2.rs_drag()

            else:
                raise io_mp.LikelihoodError(
                    "In likelihood %s. " % self.name +
                    "BAO data type %s " % self.type[i] +
                    "in %d-th line not understood" % i)

            diff.append(theo - item)


        inv_cov_ext = self.inv_cov(self.cov_data, self.error)
        chi2 = np.log(np.dot(np.dot(diff,inv_cov_ext),diff))

        # return ln(L)
        lkl = - 0.5 * chi2
        return lkl

    def Dv(self, cosmo, z): #Technically Dv*rd_fid/rd
        H = cosmo.Hubble(z)*consts.c/1000.0
        dv = ((consts.c*z*(1+z)**2*cosmo.angular_distance(z))/H)**(1/3)
        return(dv)

    def inv_cov(self, cov, errors):
        inv_cov = np.linalg.inv(cov)
        inv_cov_ext = []
        i = 0
        while i<len(errors):
            i = i + 1
            inv_cov_ext = np.append(inv_cov, 
                        np.zeros((len(errors), 
                        len(cov))), axis=0)
            
        inv_cov = np.transpose(inv_cov_ext)
        inv_cov_ext =[]
        i = 0
       
        while i<len(errors):
            i = i + 1
            inv_cov_ext = np.append(inv_cov, 
                        np.zeros((len(errors), 
                        inv_cov.shape[1])), axis=0)
            
        inv_cov_ext = np.transpose(inv_cov_ext)
        i = 1
        for item in errors:
            inv_cov_ext[-i][-i] = item**(-2)
            i = i +1 
            
        return(inv_cov_ext)

        def inv_cov_test(self, cov, errors):
        inv_cov = np.linalg.inv(cov)
        inv_cov_ext = []
        i = 0
        while i<len(errors):
            i = i + 1
            inv_cov_ext = np.append(inv_cov, 
                        np.zeros((len(errors), 
                        len(cov))), axis=0)
            
        inv_cov = np.transpose(inv_cov_ext)
        inv_cov_ext =[]
        i = 0
       
        while i<len(errors):
            i = i + 1
            inv_cov_ext = np.append(inv_cov, 
                        np.zeros((len(errors), 
                        inv_cov.shape[1])), axis=0)
            
        inv_cov_ext = np.transpose(inv_cov_ext)
        i = 1
        for item in errors:
            inv_cov_ext[-i][-i] = item**(-2)
            i = i +1 

        #Act as if there were no correlations    
        for n, row in enumerate(inv_cov_ext):
            for m, entry in enumerate(row):
                if n == m:
                    pass
                else:
                    inv_cov_ext[n][m] = 0
            
        return(inv_cov_ext)
