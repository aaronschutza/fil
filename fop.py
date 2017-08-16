#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 15:22:25 2017

@author: toffo
modified: ams7 june 2017
"""
#import cvxpy as cvx
import numpy as np
import scipy.sparse
import math
import time
import matplotlib.pyplot as plt

datafolder ='t89_series/'
dataname = 't89_kp2_tilt0.00'
neigs=10
omegas = []
xaxis = []

for nnn in range(1,100):
    datanum = str(nnn).zfill(3)
    try:
        data = np.loadtxt(datafolder+dataname+'.'+datanum+'.fop',skiprows=0)
    except:
        break
    fid = open(datafolder+dataname+'.'+datanum,"r")
    linedata = fid.readlines()
    fid.close()
    dimj = int(linedata[1])
    x=np.zeros(dimj)
    z=np.zeros(dimj)
    s=np.zeros(dimj-2)
    for i in range(4,dimj+4):
        line = [float(l) for l in linedata[i].split()]
        x[i-4] = line[0]
        z[i-4] = line[1]

    xe = x[dimj/2]

    ndim = dimj-2
    ncols = 2*ndim
    n=np.size(data)
    data.reshape((n/ncols,ncols))

    for n in range(1,dimj-2):
        s[n] = s[n-1] + np.sqrt((x[n]-x[n-1])**2+(z[n]-z[n-1])**2)
    smax=np.amax(s)
    for n in range(0,dimj-2):
        s[n] = s[n]-0.5*smax




    # #@profile(print_stats=10, dump_stats=True)
    # def sdpMethod(A):
    #     print "Testing Input Matrix..."
    #     print "\tSymmetric?\t{}".format(symQ(A));
    #     print "\tPSD? \t\t{}".format(psdQ(A));

    #     dims = A.shape;
    #     n = dims[0]
    #     X = cvx.Semidef(n)

    #     #the jumpstart below seems to take off maybe 15% ish
    #     # m = cvx.Parameter(sign="positive",value=.1)
    #     # primary = .5*cvx.sum_squares(X-A);
    #     # regularization = m*cvx.square(cvx.norm(X,2));
    #     # obj = cvx.Minimize(primary+regularization);
    #     # #constraints = [X-X.T==0,X>>0];
    #     # prob = cvx.Problem(obj);
    #     # prob.solve(solver="SCS",warm_start=False)

    #     obj = cvx.Minimize(0);
    #     constraints = [X==scipy.sparse.lil_matrix((A+A.transpose())/float(2))];
    #     #constraints = [X-X.T==0,X-A*A.transpose()==0]; different initial starting conditions
    #     prob = cvx.Problem(obj,constraints);
    #     prob.solve(solver="SCS",warm_start=False,verbose=False)
    #     print "\n"
    #     alpha = cvx.Parameter(sign="positive",value=10)
    #     beta = cvx.Parameter(sign="positive",value=.5)
    #     primary = .5*cvx.sum_squares(X-A);
    #     #2 norm forces low eigenvalues, 1 norm forces sparcity
    #     #regularization = 0;
    #     regularization = alpha*cvx.square(cvx.norm(X,2))
    #     #regularization = beta*cvx.square(cvx.norm(X,1));
    #     #regularization = alpha*cvx.square(cvx.norm(X,2))+beta*cvx.square(cvx.norm(X,1));
    #     obj = cvx.Minimize(primary+regularization);
    #     prob.objective = obj;
    #     prob.constraints = [];

    #     prob.solve(solver="SCS",warm_start=True,verbose=True,max_iters=100)

    #     Anew = scipy.sparse.coo_matrix(X.value);

    #     print "Testing New Solution..."
    #     print "\tSymmetric?\t{}".format(symQ(Anew));
    #     print "\tPSD? \t\t{}".format(psdQ(Anew));
    #     print "\tObj Val: {}".format(prob.value)
    #     return Anew;


    # def symQ(A):
    #     return (A!=A.transpose()).nnz==0

    # def psdQ(A):
    #     w = scipy.sparse.linalg.eigs(A,return_eigenvectors=False,k=1)
    #     return min(w)>=0;

    # def psdQ(A):
    #     try:
    #         scipy.linalg.cholesky(A);
    #         return True
    #     except scipy.linalg.LinAlgError as err:
    #         if 'Matrix is not positive definite' in err.message:
    #             return False



    # times = [];
    # for i in range(1,2):
    #     print "_____________________________________________________________________"
    #     n = 100*i;
    #     print "Matrix in R^({0}x{0})".format(n);
    #     A = scipy.sparse.random(n,n,density=.05,format='coo');
    #     A = 5*np.dot(A,np.transpose(A))+5*.05*n*scipy.sparse.identity(n,format='coo');

    #     E = 1/float(100)*scipy.sparse.random(n,n,density = .01,format='coo');

    #     B = scipy.sparse.coo_matrix(n,n);
    #     B = A + E;
    #     B = scipy.sparse.lil_matrix(B);

    #     start = time.clock();
    #     Anew = sdpMethod(B)

    #     end = time.clock();
    #     times.append(end-start)
    #     print "Total Time: {}\n".format(end-start);

    #     print "||Anew-A||_{0} = {1}".format("F",scipy.sparse.linalg.norm(Anew-A))
    #     print "||Anew||_{0} = {1}".format("1",scipy.sparse.linalg.norm(Anew,1))
    #     print "||Anew - Anew.T||_{0} = {1}".format("F",scipy.sparse.linalg.norm(Anew-Anew.transpose()))
    #     print "||Anew - Anew.T||_{0} = {1}".format("1",scipy.sparse.linalg.norm(Anew-Anew.transpose(),1))


    #     w = scipy.sparse.linalg.eigs(A,return_eigenvectors=False)
    #     print "Eigenvalues of A:"
    #     print "Max: {0}\t\tMin:{1}".format(max(w),min(w));

    #     w = scipy.sparse.linalg.eigs(Anew,return_eigenvectors=False)
    #     print "Eigenvalues of Anew:"
    #     print "Max: {0}\t\tMin:{1}".format(max(w),min(w));

    #     w = scipy.sparse.linalg.eigs(Anew-A,return_eigenvectors=False)
    #     print "Eigenvalues of Anew-A:"
    #     print "Max: {0}\t\tMin:{1}".format(max(w),min(w));

    # plt.scatter(range(10,1001,10),times)
    # plt.plot(range(10,1001,10), times)
    # plt.show()


    #AAA = scipy.sparse.coo_matrix(data)

    #AAA = sdpMethod(AAA)
    #AAA = 0.5*(AAA.transpose()+AAA)

    #data = AAA.toarray()

    eig_vals, eig_vecs = np.linalg.eig(data)

    eig_vals_sorted = np.sort_complex(eig_vals) 
    eig_vecs_sorted = eig_vecs[:, eig_vals.argsort()]

    
    eig = np.zeros(neigs,'complex')
    eigvp = np.zeros((neigs,ndim))
    eigvk = np.zeros((neigs,ndim))

    for n in range(neigs):
        eig[n]    = eig_vals_sorted[-1-n]
        eigvp[n,:] = eig_vecs_sorted[0:ndim,-1-n,]
        eigvk[n,:] = eig_vecs_sorted[ndim:,-1-n]
        
    omega = np.sqrt(-eig/8.52e4)
    period = 2*np.pi/omega

    omegas.append(omega.tolist())
    xaxis.append(abs(xe))
    # for i in range(neigs):
    #     titl1 = '$\Omega=$'+str(omega[i])
    #     titl2 = ' T='+str(period[i])+' sec'
    #     plt.subplot(neigs,2,2*i+1)
    #     plt.plot(s,eigvk[i,:])
    #     plt.title(titl1, fontsize=9)
    #     if i == neigs-1:
    #         plt.xlabel('s')
    #     plt.ylabel('perp')
    #     plt.grid()
    #     if i != neigs-1:
    #         plt.xticks([])
    #     plt.yticks([])
    #     plt.autoscale(enable=True, axis='x', tight=True)
    #     plt.subplot(neigs,2,2*i+2)
    #     plt.plot(s,eigvp[i,:])
    #     plt.ylabel('parl')
    #     if i == neigs-1:
    #         plt.xlabel('s')
    #     plt.title(titl2, fontsize=9)
    #     plt.grid()
    #     if i != neigs-1:
    #         plt.xticks([])
    #     plt.yticks([])
    #     plt.autoscale(enable=True, axis='x', tight=True)
    # plt.suptitle('Eigenvectors for x='+str(xe)+' Re')
    # plt.tight_layout()
    # plt.subplots_adjust(top=0.9)
    # plt.subplots_adjust(hspace=1.0)
    # plt.show()

for n in range(neigs):
    oaxis = []
    for i in range(len(xaxis)):
        oaxis.append(omegas[i][n])
    plt.plot(xaxis,oaxis)
plt.xlabel('Xe [Re]')
plt.ylabel('$\Omega$ [rad/s]')
plt.title('Eigenvalues vs. X axis')
plt.show()




