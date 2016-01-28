#!/usr/bin/env python
'''This script will profile several of the functions in the sortseq package,
sort the results by cumulative time spent on each function and then print the
results to file. The 'main' function is the script that is being targeted.
The ones that are currently targeted are profile_info using the nsb estimator,
learn_matrix using LS, learn_matrix using IM and 5 iterations, and predictiveinfo.'''

from __future__ import division
import argparse
import numpy as np
import sys
import pandas as pd
import sst.qc as qc
import sst.io as io
import os
import sst.profile_ct as profile_ct
import pdb
from sst import SortSeqError
import cProfile
import sst.profile_info as profile_info
import sst.learn_matrix as learn_matrix
import sst.predictiveinfo as predictiveinfo
import pstats

#load in data sets for the test, we will just use the sort-seq crp-wt set

df = io.load_dataset('input/mpra.txt')
model_df = io.load_model('input/mpra_model')

#Profile profile_info
#stats_fn = 'Profile_profile_info'
#stats_fn_hr = 'Profile_profile_info_hr'
#Profile.run('''profile_info.main(df,method='nsb')''',stats_fn)

#Reformat and print to human readable profile
#p = pstats.Stats(stats_fn,stream=open(stats_fn_hr,'w'))
#p.strip_dirs()
#p.sort_stats('cumtime')
#p.print_stats()

df_copy = df.copy()
#profile learn_matrix lm=LS
stats_fn = 'profile/Profile_learn_matrix_LS_mpra'
stats_fn_hr = 'profile/Profile_learn_matrix_LS_hr_mpra'
cProfile.run('''learn_matrix.main(df_copy,'dna','LS')''',stats_fn)

#Reformat and print to human readable profile
p = pstats.Stats(stats_fn,stream=open(stats_fn_hr,'w'))
p.strip_dirs()
p.sort_stats('cumtime')
p.print_stats()

'''now do the same thing to print callees data as well. This will show which functions
called each function, which could help.'''

p = pstats.Stats(stats_fn,stream=open(stats_fn_hr + '_callees','w'))
p.strip_dirs()
p.sort_stats('cumtime')
p.print_callees()

df_copy = df.copy()
#profile learn_matrix lm=IM
stats_fn = 'profile/Profile_learn_matrix_IM_mpra'
stats_fn_hr = 'profile/Profile_learn_matrix_IM_hr_mpra'
cProfile.run('''learn_matrix.main(df_copy,'dna','memsaver',iteration=5,initialize='Rand',burnin=0)''',stats_fn)

#Reformat and print to human readable profile
p = pstats.Stats(stats_fn,stream=open(stats_fn_hr,'w'))
p.strip_dirs()
p.sort_stats('cumtime')
p.print_stats()

'''now do the same thing to print callees data as well. This will show which functions
called each function, which could help.'''

p = pstats.Stats(stats_fn,stream=open(stats_fn_hr + '_callees','w'))
p.strip_dirs()
p.sort_stats('cumtime')
p.print_callees()

df_copy = df.copy()
#profile predictiveinfo
stats_fn = 'profile/Profile_predictiveinfo_mpra'
stats_fn_hr = 'profile/Profile_predictiveinfo_hr_mpra'
cProfile.run('''predictiveinfo.main(df_copy,model_df)''',stats_fn)

#Reformat and print to human readable profile
p = pstats.Stats(stats_fn,stream=open(stats_fn_hr,'w'))
p.strip_dirs()
p.sort_stats('cumtime')
p.print_stats()

'''now do the same thing to print callees data as well. This will show which functions
called each function, which could help.'''

p = pstats.Stats(stats_fn,stream=open(stats_fn_hr + '_callees','w'))
p.strip_dirs()
p.sort_stats('cumtime')
p.print_callees()
