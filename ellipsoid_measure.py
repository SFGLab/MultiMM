from scipy.spatial import ConvexHull, convex_hull_plot_2d
from numpy.linalg import eig, inv
import numpy as np
import pandas as pd
import os
import math 
import sys
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from numpy.linalg import eig, inv
from points_io import save_points_as_pdb, point_reader

def ls_ellipsoid(xx,yy,zz):                                  
    #finds best fit ellipsoid. Found at http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html
    #least squares fit to a 3D-ellipsoid
    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz  = 1
    #
    # Note that sometimes it is expressed as a solution to
    #  Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz  = 1
    # where the last six terms have a factor of 2 in them
    # This is in anticipation of forming a matrix with the polynomial coefficients.
    # Those terms with factors of 2 are all off diagonal elements.  These contribute
    # two terms when multiplied out (symmetric) so would need to be divided by two
    
    # change xx from vector of length N to Nx1 matrix so we can use hstack
    x = xx[:,np.newaxis]
    y = yy[:,np.newaxis]
    z = zz[:,np.newaxis]
    
    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz = 1
    J = np.hstack((x*x,y*y,z*z,x*y,x*z,y*z, x, y, z))
    K = np.ones_like(x) #column of ones
    
    #np.hstack performs a loop over all samples and creates
    #a row in J for each x,y,z sample:
    # J[ix,0] = x[ix]*x[ix]
    # J[ix,1] = y[ix]*y[ix]
    # etc.
    
    JT=J.transpose()
    JTJ = np.dot(JT,J)
    InvJTJ=np.linalg.inv(JTJ);
    ABC= np.dot(InvJTJ, np.dot(JT,K))

    # Rearrange, move the 1 to the other side
    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz - 1 = 0
    #    or
    #  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz + J = 0
    #  where J = -1
    eansa=np.append(ABC,-1)

    return (eansa)

def polyToParams3D(vec,printMe):                             
    #gets 3D parameters of an ellipsoid. Found at http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html
    # convert the polynomial form of the 3D-ellipsoid to parameters
    # center, axes, and transformation matrix
    # vec is the vector whose elements are the polynomial
    # coefficients A..J
    # returns (center, axes, rotation matrix)
    
    #Algebraic form: X.T * Amat * X --> polynomial form
    
    if printMe: print('\npolynomial\n',vec)
    
    Amat=np.array(
    [
    [ vec[0],     vec[3]/2.0, vec[4]/2.0, vec[6]/2.0 ],
    [ vec[3]/2.0, vec[1],     vec[5]/2.0, vec[7]/2.0 ],
    [ vec[4]/2.0, vec[5]/2.0, vec[2],     vec[8]/2.0 ],
    [ vec[6]/2.0, vec[7]/2.0, vec[8]/2.0, vec[9]     ]
    ])
    
    if printMe: print('\nAlgebraic form of polynomial\n',Amat)
    
    #See B.Bartoni, Preprint SMU-HEP-10-14 Multi-dimensional Ellipsoidal Fitting
    # equation 20 for the following method for finding the center
    A3=Amat[0:3,0:3]
    A3inv=inv(A3)
    ofs=vec[6:9]/2.0
    center=-np.dot(A3inv,ofs)
    if printMe: print('\nCenter at:',center)
    
    # Center the ellipsoid at the origin
    Tofs=np.eye(4)
    Tofs[3,0:3]=center
    R = np.dot(Tofs,np.dot(Amat,Tofs.T))
    if printMe: print('\nAlgebraic form translated to center\n',R,'\n')
    
    R3=R[0:3,0:3]
    R3test=R3/R3[0,0]
    # print('normed \n',R3test)
    s1=-R[3, 3]
    R3S=R3/s1
    (el,ec)=eig(R3S)
    
    recip=1.0/np.abs(el)
    axes=np.sqrt(recip)
    if printMe: print('\nAxes are\n',axes  ,'\n')
    
    inve=inv(ec) #inverse is actually the transpose here
    if printMe: print('\nRotation matrix\n',inve)
    return (center,axes,inve)


def calc_ellipsoid_ratio(V, result="ratio"):
    """
    V: np.array
    result: either "ratio" or "volume"
    """
    x = V[:,0]
    y = V[:,1]
    z = V[:,2]

    #get convex hull
    surface  = np.stack((x,y,z), axis=-1)
    hullV    = ConvexHull(surface)
    lH       = len(hullV.vertices)
    hull     = np.zeros((lH,3))
    for i in range(len(hullV.vertices)):
        hull[i] = surface[hullV.vertices[i]]
    hull     = np.transpose(hull)         

    #fit ellipsoid on convex hull
    eansa            = ls_ellipsoid(hull[0],hull[1],hull[2]) #get ellipsoid polynomial coefficients
    center,axes,inve = polyToParams3D(eansa,False)   #get ellipsoid 3D parameters
    
    if result == "ratio":
        return np.max(axes)/np.min(axes)
    elif result == "volume":
        return 4/3*np.pi*axes[0]*axes[1]*axes[2]
    else:
        raise(Exception("Wrong result provided: {}. Expected: ratio or volume.".format(result)))
    
    
def calc_ellipsoid_ratio2(V, n_e = 20, o_e = 0.2, result="ratio"):
    """
    V: np.array
    result: either "ratio" or "volume"
    n_e: number of small ellipsoids to compute ellipsoid measure
    o_e: overlap between consequtive small ellipsoids. Should be between 0 and 1
    """
    n = V.shape[0]
    
    s_e = n/(n_e*(1-o_e)+o_e)
    beginnings = [i*(1-o_e)*s_e for i in range(n_e)]
    
    ratios = []
    for b in beginnings:
        ratios.append(calc_ellipsoid_ratio(V[int(b):int(b+s_e),:], result=result))
    
    return np.mean(ratios) 







