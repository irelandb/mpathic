#!/usr/bin/env python

from __future__ import division
import os
#Our standard Modules
import argparse
import numpy as np
import scipy as sp
import sys
#Our miscellaneous functions
import pandas as pd
import sst.utils as utils
from sklearn import linear_model
'''commands to run all analysis. WARNING, this will take ~ 30 hours total.
Commands should probably be run individually'''

'''MCMC analysis'''
#learn CRP models
os.system('''sst learn_matrix -lm IM -i crp-wt_formatted.txt -o crp-wt_MCMC_crp -s 3 -e 25''')
os.system('''sst learn_matrix -lm LS -i crp-wt_formatted.txt -o crp-wt_LS_crp -s 3 -e 25''')
os.system('''sst learn_matrix -lm ER -i crp-wt_formatted.txt -o crp-wt_ER_crp -s 3 -e 25''')
os.system('''sst learn_matrix -lm IM -i full-wt_formatted.txt -o full-wt_MCMC_crp -s 3 -e 25''')
os.system('''sst learn_matrix -lm LS -i full-wt_formatted.txt -o full-wt_LS_crp -s 3 -e 25''')
os.system('''sst learn_matrix -lm ER -i full-wt_formatted.txt -o full-wt_ER_crp -s 3 -e 25''')
os.system('''sst learn_matrix -lm IM -i full-500_formatted.txt -o full-500_MCMC_crp -s 3 -e 25''')
os.system('''sst learn_matrix -lm LS -i full-500_formatted.txt -o full-500_LS_crp -s 3 -e 25''')
os.system('''sst learn_matrix -lm ER -i full-500_formatted.txt -o full-500_ER_crp -s 3 -e 25''')
os.system('''sst learn_matrix -lm IM -i full-150_formatted.txt -o full-150_MCMC_crp -s 3 -e 25''')
os.system('''sst learn_matrix -lm LS -i full-150_formatted.txt -o full-150_LS_crp -s 3 -e 25''')
os.system('''sst learn_matrix -lm ER -i full-150_formatted.txt -o full-150_ER_crp -s 3 -e 25''')
os.system('''sst learn_matrix -lm IM -i full-0_formatted.txt -o full-0_MCMC_crp -s 3 -e 25''')
os.system('''sst learn_matrix -lm LS -i full-0_formatted.txt -o full-0_LS_crp -s 3 -e 25''')
os.system('''sst learn_matrix -lm ER -i full-0_formatted.txt -o full-0_ER_crp -s 3 -e 25''')

#learn CRP neighbor models
os.system('''sst learn_matrix -lm IM -i crp-wt_formatted.txt -o crp-wt_MCMC_crp_NBR -s 3 -e 25 -mt NBR''')
os.system('''sst learn_matrix -lm IM -i full-wt_formatted.txt -o full-wt_MCMC_crp_NBR -s 3 -e 25 -mt NBR''')
os.system('''sst learn_matrix -lm IM -i full-500_formatted.txt -o full-500_MCMC_crp_NBR -s 3 -e 25 -mt NBR''')
os.system('''sst learn_matrix -lm IM -i full-150_formatted.txt -o full-150_MCMC_crp_NBR -s 3 -e 25 -mt NBR''')
os.system('''sst learn_matrix -lm IM -i full-0_formatted.txt -o full-0_MCMC_crp_NBR -s 3 -e 25 -mt NBR''')

#learn RNAP models
os.system('''sst learn_matrix -lm IM -i rnap-wt_formatted.txt -o rnap-wt_MCMC_rnap -s 37 -e 71''')
os.system('''sst learn_matrix -lm LS -i rnap-wt_formatted.txt -o rnap-wt_LS_rnap -s 37 -e 71''')
os.system('''sst learn_matrix -lm ER -i rnap-wt_formatted.txt -o rnap-wt_ER_rnap -s 37 -e 71''')
os.system('''sst learn_matrix -lm IM -i full-wt_formatted.txt -o full-wt_MCMC_rnap -s 37 -e 71''')
os.system('''sst learn_matrix -lm LS -i full-wt_formatted.txt -o full-wt_LS_rnap -s 37 -e 71''')
os.system('''sst learn_matrix -lm ER -i full-wt_formatted.txt -o full-wt_ER_rnap -s 37 -e 71''')
os.system('''sst learn_matrix -lm IM -i full-500_formatted.txt -o full-500_MCMC_rnap -s 37 -e 71''')
os.system('''sst learn_matrix -lm LS -i full-500_formatted.txt -o full-500_LS_rnap -s 37 -e 71''')
os.system('''sst learn_matrix -lm ER -i full-500_formatted.txt -o full-500_ER_rnap -s 37 -e 71''')
os.system('''sst learn_matrix -lm IM -i full-150_formatted.txt -o full-150_MCMC_rnap -s 37 -e 71''')
os.system('''sst learn_matrix -lm LS -i full-150_formatted.txt -o full-150_LS_rnap -s 37 -e 71''')
os.system('''sst learn_matrix -lm ER -i full-150_formatted.txt -o full-150_ER_rnap -s 37 -e 71''')
os.system('''sst learn_matrix -lm IM -i full-0_formatted.txt -o full-0_MCMC_rnap -s 37 -e 71''')
os.system('''sst learn_matrix -lm LS -i full-0_formatted.txt -o full-0_LS_rnap -s 37 -e 71''')
os.system('''sst learn_matrix -lm ER -i full-0_formatted.txt -o full-0_ER_rnap -s 37 -e 71''')

#learn RNAP neighbor models
os.system('''sst learn_matrix -lm IM -i crp-wt_formatted.txt -o crp-wt_MCMC_crp_NBR -s 37 -e 71 -mt NBR''')
os.system('''sst learn_matrix -lm IM -i full-wt_formatted.txt -o full-wt_MCMC_crp_NBR -s 37 -e 71 -mt NBR''')
os.system('''sst learn_matrix -lm IM -i full-500_formatted.txt -o full-500_MCMC_crp_NBR -s 37 -e 71 -mt NBR''')
os.system('''sst learn_matrix -lm IM -i full-150_formatted.txt -o full-150_MCMC_crp_NBR -s 37 -e 71 -mt NBR''')
os.system('''sst learn_matrix -lm IM -i full-0_formatted.txt -o full-0_MCMC_crp_NBR -s 37 -e 71 -mt NBR''')

#learn DMS models of the two data sets
os.system('''sst learn_matrix -lm IM -i dms_1_formatted -o dms_1_MCMC''')
os.system('''sst learn_matrix -lm LS -i dms_1_formatted -o dms_1_LS''')
os.system('''sst learn_matrix -lm ER -i dms_1_formatted -o dms_1_ER''')
os.system('''sst learn_matrix -lm IM -i dms_2_formatted -o dms_2_MCMC''')
os.system('''sst learn_matrix -lm LS -i dms_2_formatted -o dms_2_LS''')
os.system('''sst learn_matrix -lm ER -i dms_2_formatted -o dms_2_ER''')

#learn mpra models 
os.system('''sst learn_matrix -lm IM -i CRE_100uM_ds_formatted -o CRE_100uM_ds_MCMC''')
os.system('''sst learn_matrix -lm LS -i CRE_100uM_ds_formatted -o CRE_100uM_ds_LS''')
os.system('''sst learn_matrix -lm ER -i CRE_100uM_ds_formatted -o CRE_100uM_ds_ER''')
os.system('''sst learn_matrix -lm IM -i CRE_100uM_test_formatted -o CRE_100uM_test_MCMC''')
os.system('''sst learn_matrix -lm LS -i CRE_100uM_test_formatted -o CRE_100uM_test_LS''')
os.system('''sst learn_matrix -lm ER -i CRE_100uM_test_formatted -o CRE_100uM_test_ER''')

#learn mpra Neighbor models
os.system('''sst learn_matrix -lm IM -i CRE_100uM_ds_formatted -o CRE_100uM_ds_MCMC_NBR -mt NBR''')
os.system('''sst learn_matrix -lm IM -i CRE_100uM_test_formatted -o CRE_100uM_test_MCMC_NBR -mt NBR''')
'''calculate the predictive information'''

#CRP models
#run all predictive info calculations
models = ['crp-wt_MCMC_crp','full-wt_MCMC_crp','full-500_MCMC_crp','full-150_MCMC_crp','full-0_MCMC_crp']
modeltype = ['_MCMC','_IM','_LS']
datasets = ['crp-wt_formatted.txt','full-wt_formatted.txt','full-500_formatted.txt','full-150_formatted.txt','full-0_formatted.txt']
MI = np.zeros([5,5])
for i in models:
    for mt in modeltype:
        for z in datasets:
            outname = 'info_' + i + z 
            os.system('''sst predictiveinfo -m %s -ds %s -o %s''' % (i + mt,z,'MI_' + i + mt + '_' + z))
        
        

