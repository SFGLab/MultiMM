import pandas as pd
import numpy as np
import torch

def amplify(structure, scale=10):
    return [(s[0]*scale, s[1]*scale, s[2]*scale) for s in structure]

def get_orig_helix(n=40, r=3.25, a=0.2):
    points = []
    for t in np.linspace(0, 4*np.pi, n):
        points.append((r*np.cos(t)-r, r*np.sin(t), a*t))
    return points

def orthonormalize(vectors):
    """
        Orthonormalizes the vectors using gram schmidt procedure.

        Parameters:
            vectors: torch tensor, size (dimension, n_vectors)
                    they must be linearly independant
        Returns:
            orthonormalized_vectors: torch tensor, size (dimension, n_vectors)
    """
    assert (vectors.size(1) <= vectors.size(0)), 'number of vectors must be smaller or equal to the dimension'
    orthonormalized_vectors = torch.zeros_like(vectors)
    orthonormalized_vectors[:, 0] = vectors[:, 0] / torch.norm(vectors[:, 0], p=2)

    for i in range(1, orthonormalized_vectors.size(1)):
        vector = vectors[:, i]
        V = orthonormalized_vectors[:, :i]
        PV_vector= torch.mv(V, torch.mv(V.t(), vector))
        orthonormalized_vectors[:, i] = (vector - PV_vector) / torch.norm(vector - PV_vector, p=2)

    return orthonormalized_vectors

def move_structure_to(struct, p1, p2, x0=np.array([None])):
    if p1[0]==p2[0] and p1[1]==p2[1] and p1[2]==p2[2]:
        raise(Exception("Starting point and the ending points must be different!"))
    v_x = [1,0,0]
    v_y = [0,1,0]
    v_z = [0,0,1]
    w_z = [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]
    if w_z[0] != 0:
        w_x = [-w_z[1]/w_z[0], 1, 0]
        w_y = [-w_z[2]/w_z[0], 0, 1]
    elif w_z[1] != 0:
        w_x = [1, -w_z[0]/w_z[1], 0]
        w_y = [0, -w_z[2]/w_z[1], 1]
    else:
        w_x = [1, 0, -w_z[0]/w_z[2]]
        w_y = [0, 1, -w_z[1]/w_z[2]]
    A = torch.transpose(torch.tensor([w_z, w_x, w_y]),0,1)
    A_norm = np.array(orthonormalize(A))
    w_x = list(A_norm[:,1])
    w_y = list(A_norm[:,2])
    w_z = list(A_norm[:,0])
    if x0.all()==None: x0=p1
    new_helix = []
    for p in struct:
        new_helix.append((x0[0]+p[0]*w_x[0]+p[1]*w_y[0]+p[2]*w_z[0],
                         x0[1]+p[0]*w_x[1]+p[1]*w_y[1]+p[2]*w_z[1],
                         x0[2]+p[0]*w_x[2]+p[1]*w_y[2]+p[2]*w_z[2]))
    return new_helix

def rot_mat(th,axis='z'):
    if axis=='x':
        mat = np.array([[1,0,0],
                        [0,np.cos(th),-np.sin(th)],
                        [0,np.sin(th),np.cos(th)]])
    elif axis=='y':
        mat = np.array([[np.cos(th),0,np.sin(th)],
                        [0,1,0],
                        [-np.sin(th),0,np.cos(th)]])
    elif axis=='z':
        mat = np.array([[np.cos(th),-np.sin(th),0],
                        [np.sin(th), np.cos(th),0],
                        [0,0,1]])
    return mat

def add_helixes(points,entry,exit):
    helix = get_orig_helix(n=40, r=3.25, a=0.2)
    nuc_pos = list()
    for en, ex in zip(entry,exit):
        p1, p2 = points[en], points[ex]
        new_helix = move_structure_to(helix, p1, p2)
        points[en:ex] = new_helix
        nuc_pos.append(np.average(new_helix,axis=0))
    nuc_pos = np.array(nuc_pos)
    return points, nuc_pos