#!/usr/bin/env python
import time
import mpathic.simulate_library
import mpathic.fast as fast
import mpathic.qc as qc
from mpathic.profile_mut import main as profile_mut
from mpathic.simulate_library import main as simulate_library
import numpy as np
from scipy.sparse import csr, csr_matrix, lil_matrix
from mpathic import SortSeqError
import pdb
import sys
import pandas as pd
import mpathic.numerics as numerics
from numpy.random import randn
import mpathic.Models as Models
import time

# Create sequences to test this on
wtseq = 'AAAAAAAGTGAGATGGCAATCTAATTCGGCACCCCAGGTTTTACACTTTATGCTTCCGGCTCGTATGTTGTGTGG'
L = len(wtseq)
modeltypes = ['MAT','NBR']
seqtypes = ['dna','protein']
numseqs_dict = {'dna':10000,'protein':1000}
for seqtype in seqtypes:
    for modeltype in modeltypes:
        for mutrate in [0.01,0.1,1]:
            numseqs = numseqs_dict[seqtype]
            dataset_df = simulate_library(wtseq,numseq=numseqs,mutrate=mutrate,tags=True,\
                dicttype=seqtype)
            seqarray = numerics.dataset2seqarray(dataset_df,\
                modeltype=modeltype)
            mutarray, wtrow = numerics.dataset2mutarray(dataset_df,\
                modeltype=modeltype)

            # Print compression results
            seqarray_size = numerics.nbytes(seqarray)
            mutarray_size = numerics.nbytes(mutarray)

            # Create matrix for random model
            alphabet = qc.seqtype_to_alphabet_dict[seqtype]
            C = len(alphabet)
            num_rows = (L-1) if modeltype=='NBR' else L
            num_cols = C**2 if modeltype=='NBR' else C
            modelmatrix = randn(num_rows,num_cols)

            # Create model dataframe
            val_cols = qc.model_parameters_dict[(modeltype,seqtype)]
            model_df = pd.DataFrame(modelmatrix,columns=val_cols)
            model_df['pos'] = range(num_rows)
            model_df = qc.validate_model(model_df)

            # Evaluate model on sequences the standard way
            if modeltype == 'MAT':
                model_obj = Models.LinearModel(model_df)
            elif modeltype == 'NBR':
                model_obj = Models.NeighborModel(model_df)
            else:
                raise SortSeqError('Unrecognized model type %s'%modeltype)

            # Evaluate model on sequences using slow method
            t = time.time()
            vals_slow = model_obj.evaluate_on_seqarray(seqarray)
            t_slow = time.time()-t

            # Evaluate model on sequences using fast method
            t = time.time()
            vals_fast = model_obj.evaluate_on_mutarray(mutarray,wtrow)
            t_fast = time.time()-t
            print '%s, %s, %.2f: compression ratio = %.1f, speedup ratio = %.1f'%\
                (modeltype, seqtype, mutrate,\
                1.*seqarray_size/mutarray_size, t_slow/t_fast)
            assert np.allclose(vals_slow, vals_fast)



