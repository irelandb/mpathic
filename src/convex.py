import numpy as np
from cvxopt import solvers
import cvxopt
import scipy.sparse
import warnings

def ridge_objective(params, mean_bins, Xmat, alpha):
    '''
    Computes ridge objective function (sum of squares of error with L2
    penalty.
    '''
    predicted = np.dot(Xmat, params)
    r_obj = ((y - mean_bins)**2).sum()
    r_obj += alpha * (params**2).sum()

    return r_obj

def test_iter(A, B):
    m,n1 = A.shape
    n2 = B.shape[1]
    Cshape = (m, n1*n2)
    with warnings.catch_warnings():
        #ignore depreciation warnings
        warnings.simplefilter("ignore")
    
        #initialize our output sparse matrix data objects
        data = np.empty((m,),dtype=object)
        col =  np.empty((m,),dtype=object)
        row =  np.empty((m,),dtype=object)
        #do multiplication for each row
        for i,(a,b) in enumerate(zip(A, B)):
            #column indexes
            col1 = a.indices * n2
            col[i] = (col1[:,None]+b.indices).flatten()
            row[i] = np.full((a.nnz*b.nnz,), i)
            #all data will be 1's, as it is only true or false data
            data[i] = np.ones(len(col[i]))
    data = np.concatenate(data)
    col = np.concatenate(col)
    row = np.concatenate(row)
    return scipy.sparse.coo_matrix((data,(row,col)),shape=Cshape)

def ridge_gradient(params, mean_bins, Xmat, alpha):
    '''
    The gradiant of the ridge regression objective.
    '''
    grad = 2*(np.dot(np.dot(Xmat.T, Xmat), params) - np.dot(Xmat.T, y))
    grad += 2*(alpha * params)
    return grad

def ridge_hess(params, mean_bins, Xmat, alpha):
    '''
    The hessian of the ridge regression objective.
    '''
    hessian = np.dot(Xmat.T, Xmat) + np.eye(2) * alpha
    return hessian

def test_iter(A, B):
    m,n1 = A.shape
    n2 = B.shape[1]
    Cshape = (m, n1*n2)
    with warnings.catch_warnings():
        #ignore depreciation warnings
        warnings.simplefilter("ignore")
    
        #initialize our output sparse matrix data objects
        data = np.empty((m,),dtype=object)
        col =  np.empty((m,),dtype=object)
        row =  np.empty((m,),dtype=object)
        #do multiplication for each row
        for i,(a,b) in enumerate(zip(A, B)):
            #column indexes
            col1 = a.indices * n2
            col[i] = (col1[:,None]+b.indices).flatten()
            row[i] = np.full((a.nnz*b.nnz,), i)
            #all data will be 1's, as it is only true or false data
            data[i] = np.ones(len(col[i]))
    data = np.concatenate(data)
    col = np.concatenate(col)
    row = np.concatenate(row)
    return scipy.sparse.coo_matrix((data,(row,col)),shape=Cshape)

def convex_opt_agorithm(s,N0, Nsm,tm,alpha=1):
    '''this function will use cvxopt to perform convex optimization
        on a poissonian objective for the PR learning method.'''
    bins = tm.shape[1]
    N0_matrix = np.matrix(N0)
    tm_matrix = np.matrix(tm)
    Nsm_matrix = np.matrix(Nsm)
    tm_matrix_squared = np.matrix(np.multiply(tm,tm))
    i,c = s.shape
    #create s matrix, the elements of this matrix are delta_s@i * delta s2@j
    #we will need this for hessian.
    s_hessian_mat = scipy.sparse.csr_matrix(test_iter(s,s))
    def F(x=None, z=None):
        if x is None: return 0, cvxopt.matrix(0.0, (c+bins,1))
        gm = np.transpose(np.array(x[c:]))
        y = s*x[:c]
        w = np.add(gm,y*tm_matrix)
        term2 = np.array(N0)*np.array(np.exp(w))
        term1 = np.multiply(np.array(Nsm),np.array(w))
        f = cvxopt.matrix(-np.sum(term1-term2)+alpha/2*np.sum(cvxopt.mul(x,x)))
        Df = np.zeros((1,c+bins))
        
        Nm = np.sum(Nsm,axis=0)
        
        
        Df_gm = Nm - sum(term2)
        
        Df_theta =  np.transpose((Nsm - term2)*np.transpose(tm))*s
    
        Df[0,:c] = Df_theta
        
        Df[0,c:] = Df_gm
        Df = Df - np.transpose(alpha*x)
        Df_cvx = cvxopt.matrix(-Df)
        if z is None: return f, Df_cvx
        H = np.zeros((c+bins,c+bins))
        
        #first do the theta_theta terms
                
        #first do N_0*tm**2, this will produce a ixm matrix where each term is Ni0 * tm**2.
        Inner_term = N0*tm_matrix_squared
        #now multiply by a column vector form of w, this will do the multiplication by y, and also sum over m.
        #you now have an ix1 matrix
        Inner_term2 = np.sum(np.array(Inner_term)*np.array(np.exp(w)),axis=1)
        
        #multiply by hessian matrix and sum over sequences...
     
        H_thetas_temp = Inner_term2*s_hessian_mat
        
        #now we need to reshape such to fill in these values
        H[:c,:c] = H_thetas_temp.reshape((c,c),order='F')
        
        
        
        
        #now do mixed terms with partials of gms and thetas
        #first multiply the N0 counts by the exponential term. This will give you a matrix of ixm
        Inner_term = np.array(N0)*np.array(np.matrix(np.exp(w)))
        #now use element wise multiplication to multiply by the tms, this will broadcasts their value down the rows
        Inner_term2 = np.array(tm)*Inner_term
        #now convert back to matrix form and multiply by s. This sums over sequence while multiplying by si, and
        #yeilds a matrix of dimension mxc
        H_gm_theta_temp = np.transpose(np.matrix(Inner_term2))*s
        H[c:c+bins,:c] = H_gm_theta_temp
        
        #now we are going to do the last bit, the partials with respect to the gms
        
        #we can use the same 'Inner term' as above N0*exp(gm+tm*sum(theta*s))
        #sum over s
        Temp_term = np.sum(Inner_term,axis=0)
        ident = np.identity(bins)
        #use array multiplication to broadcast multiply such that all off diagonal terms are 0 and
        #the on diagonal terms are equal to sum(Ns0*exp(gm+sum(theta*s)))
        H_gms = np.array(Temp_term)*ident
        H[c:c+bins,c:c+bins] = H_gms
        
        
        #now insure that H is symmetric across diagonal
        for q in range(c+bins):
            for k in range(q):
                H[k,q] = H[q,k]
        penalty_term = alpha*np.identity(c+bins)
        H = cvxopt.matrix(z[0]*(H+penalty_term))
                  
        
        return f, Df_cvx, H
    return solvers.cp(F)['x']

def convex_opt_agorithm_LS(s,Nsm,tm,alpha=1):
    '''This function will use cvxopt to perform convex optimization
        on the least squares objective'''
    bins = tm.shape[1]
    N0_matrix = np.matrix(N0)
    tm_matrix = np.matrix(tm)
    #get weighted mean of bin number and counts
    bin_numbers = np.array([int(q.split('_')[1]) for q in Nsm.columns])
    temp_Nsm = np.array(Nsm)*bin_numbers
    mean_bins = np.mean(temp_Nsm,axis=1)
    i,c = s.shape
    #create s matrix, the elements of this matrix are delta_s@i * delta s2@j
    #we will need this for hessian.
    s_hessian_mat = scipy.sparse.csr_matrix(test_iter(s,s))
    def F(x=None, z=None):
        if x is None: return 0, matrix(0.0, (c,1))
        y = s*x - mean_counts
        f = cvxopt.matrix(ridge_objective(x,mean_bins,s,alpha))
        Df = cvxopt.matrix(ridge_gradient(x,mean_bins,s,alpha))
        if z is None: return f, Df
        H = cvxopt.matrix(ridge_hess(x,mean_bins,s,alpha))
                  
        
        return f, Df_cvx, H
    return solvers.cp(F)['x']
