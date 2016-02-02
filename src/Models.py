#!/usr/bin/env python
'''Class container for Model Handling in the sortseq package. We currently
    only support linear and nearest neighbor models for use within the package.
    In the future, support for additional model types will be added'''
import importlib
import sys
import numpy as np
import sst.utils as utils	
import scipy as sp
import pandas as pd
import pdb
import sst.qc as qc
import sst.io as io
import time
from sst import SortSeqError


class ExpModel:
    #All exp models return log expression
    pass


class LinearModel(ExpModel):
    ''' This model generates energies of binding 
        through a linear energy matrix model'''

    def __init__(self,model_df):
        """
        Constructor takes model parameters in the form of a model dataframe
        """
        model_df = qc.validate_model(model_df.copy(),fix=True)
        seqtype, modeltype = qc.get_model_type(model_df)
        if not modeltype=='MAT':
            raise SortSeqError('Invalid modeltype: %s'%modeltype)

        seq_dict,inv_dict = utils.choose_dict(seqtype,modeltype=modeltype)
        self.seq_dict = seq_dict
        self.inv_dict = inv_dict
        self.df = model_df
        self.length = model_df.shape[0]

        # Extract matrix part of model dataframe
        headers = qc.get_cols_from_df(model_df,'vals')
        self.matrix = np.transpose(np.array(model_df[headers]))

    def genexp(self,seqs):
        return self.evaluate(seqs)

    def evaluate(self,seqs):
        # Check seqs container
        if not (isinstance(seqs,list) or isinstance(seqs,pd.Series)):
            raise SortSeqError('Sequences must be input as a list')
        # Check length
        if len(seqs[0]) != self.matrix.shape[1]:
            raise SortSeqError(\
                'Energy Matrix Length does not equal Sequence Length')

        # Compute seqmats
        t0 = time.time()
        seqmats = [utils.seq2mat(s,self.seq_dict) for s in seqs]
        t1 = time.time()

        # Compute and return values
        vals = np.array([np.sum(self.matrix*s) for s in seqmats])
        t2 = time.time()

        #print 't1-t0 = %.4f, t1-t2 = %.4f'%(t1-t0,t2-t1)
        return vals 

class NeighborModel(ExpModel):
    ''' This model generates energies of binding 
        through a nearest neighbor matrix model'''

    def __init__(self,model_df):
        """
        Constructor takes model parameters in the form of a model dataframe
        """
        model_df = qc.validate_model(model_df.copy(),fix=True)
        seqtype, modeltype = qc.get_model_type(model_df)
        if not modeltype=='NBR':
            raise SortSeqError('Invalid modeltype: %s'%modeltype)

        seq_dict,inv_dict = utils.choose_dict(seqtype,modeltype=modeltype)
        self.seq_dict = seq_dict
        self.inv_dict = inv_dict
        self.df = model_df
        self.length = model_df.shape[0]+1

        # Extract matrix part of model dataframe
        headers = qc.get_cols_from_df(model_df,'vals')
        self.matrix = np.transpose(np.array(model_df[headers]))

    def genexp(self,seqs):
        return self.evaluate(seqs)

    def evaluate(self,seqs):
        # Check seqs container
        if not (isinstance(seqs,list) or isinstance(seqs,pd.Series)):
            raise SortSeqError('Sequences must be input as a list')
        # Check length
        if len(seqs[0]) != self.matrix.shape[1]+1:
            raise SortSeqError(\
                'Energy Matrix Length does not equal Sequence Length')

        # Compute seqmats
        t0 = time.time()
        seqmats = [utils.seq2matpair(s,self.seq_dict) for s in seqs]
        t1 = time.time()

        # Compute and return values
        vals = np.array([np.sum(self.matrix*s) for s in seqmats])
        t2 = time.time()

        #print 't1-t0 = %.4f, t1-t2 = %.4f'%(t1-t0,t2-t1)

        return vals 

class RaveledModel(ExpModel):
    '''This model type generates an energy matrix model, 
        but accepts flattened matrices (for sklearn based fitting). It will not
        generate models from files or dataframes.'''

    def __init__(self,param):
          self.matrix = param
     
    def genexp(self,sparse_seqs_mat,sample_weight=None):
        if sp.sparse.issparse(sparse_seqs_mat):
            '''If modeltype is an energy matrix for repression or activation,
                this will calculate the binding energy of a sequence, which will 
                be monotonically correlated with expression.'''
            n_seqs = sparse_seqs_mat.shape[0]
            #energies = sp.zeros(n_seqs) 
            energies = np.zeros(n_seqs)
            energiestemp = sparse_seqs_mat.multiply(self.matrix).sum(axis=1)
            for i in range(0,n_seqs):
                 energies[i] = energiestemp[i]
        else:
            raise ValueError('Enter sparse sequence matrix.')
        if sample_weight:
            t_exp = np.zeros(np.sum(sample_weights))
            counter=0
            for i, sw in enumerate(sample_weight):
                t_exp[counter:counter+sample_weight[i]] = -energies[i]
                counter = counter + sample_weight[i]
            return t_exp
        else:
            return energies

class NoiseModel():

    def gennoisyexp(self,logexp):
        '''If no particular noise function is declared, then just assume no
            noise'''
        return logexp

    def genlist(self,df):
        nlist = []
        T_counts = np.sum(df['ct'])
        nexp = np.zeros(T_counts)
        counter = 0
        
        for i,seq in enumerate(df['seq']):
            exp = df['val'][i]
            counts = df['ct'][i]
            explist = [exp for z in range(counts)]
            noisyexp = self.gennoisyexp(explist)
            nlist.append(noisyexp)
            nexp[counter:counter+counts]=noisyexp
            counter = counter + counts
        return nexp,nlist
    
        

class LogNormalNoise(NoiseModel):
    '''Noise model that adds autoflourescence and then draws the noisy 
        measurement from a log normal distribution with this mean'''

    def __init__(self,npar):
        self.auto = float(npar[0])
        self.scale = float(npar[1])

    def gennoisyexp(self,logexp):
        exp = np.exp(logexp) + self.auto
        nexp = np.random.lognormal(np.log(exp),self.scale)
        return np.log(nexp)

class NormalNoise(NoiseModel):
    '''Add Gaussian Noise'''

    def __init__(self,npar):
        try:
            self.scale = float(npar[0])
        except ValueError:
            raise IOError('your input parameter must be a float')
        #Check that scale is in the correct range
        if self.scale < 0:
            raise IOError('''your input scale for normal noise must be greater
                than zero''')       

    def gennoisyexp(self,logexp):
        #set scale of normal distribution and add 1e-6 to make sure its not zero
        s = (np.abs(logexp).mean())*self.scale + 1e-6
        #add noise to each entry
        nexp = logexp + np.random.normal(
            scale=s,size=len(logexp))
        return nexp

class PoissonNoise(NoiseModel):
    '''Add Noise for mpra experiment Expression Measurements'''

    def gennoisyexp(self,df,T_LibCounts,T_mRNACounts):
        exp = np.exp(df['val'])
        libcounts = df['ct']
        weights = exp*libcounts
        meanexp = T_mRNACounts*weights/np.sum(weights)
        meanlibcounts = libcounts/np.sum(libcounts)*T_LibCounts
        noisyexp = np.random.poisson(lam=meanexp)
        noisylibcounts = np.random.poisson(lam=meanlibcounts)
        return noisylibcounts,noisyexp
