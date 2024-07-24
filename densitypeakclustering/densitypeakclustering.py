#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Script custom implements the fast search and clustering by density peaks
algorithm

Rodriguez, Alex, and Alessandro Laio. “Clustering by Fast Search and Find of Density Peaks.” Science 344, no. 6191 (June 27, 2014): 1492-96.

Created on Sunday 1 Aug 2021

@author: pgoltstein
"""



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Imports
import numpy as np
import sklearn.metrics
import matplotlib.pyplot as plt



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions

def find_clusters(X, fraction=0.05, rho_min=0.2, delta_min=0.2, weights=None, rho_x_delta_min=None, show_rho_vs_delta=False, quiet=False):
    """
    Runs the cluster detection from start to end.

    Parameters
    ----------
    X : ndarray
        Data matrix X (samples x features)
    fraction : float, optional
        Fraction of data to include for distance calculation, by default 0.05
    rho_min : float, optional
        Minimum rho for cluster to be included, by default 0.2
    delta_min : float, optional
        Minimum delta for cluster to be included, by default 0.2
    weights : ndarray, optional
        Value to scale the density contribution of each datapoint, by default None
    rho_x_delta_min : float, optional
        Minimum product of rho and delta to be included as cluster, if value given, it supersedes rho_min and delta_min, by default None
    show_rho_vs_delta : bool, optional
        Whether to show scatter plot of rho vs delta, by default False
    quiet : bool, optional
        Suppress print statements, by default False

    Returns
    -------
    list
        List that holds dictionaries with cluster stats
        [cluster = {'X':#, 'Y':#, 'rho':#, 'delta':#}, ...]
    """

    # Output list
    clusters = []

    # Cluster detection cascade
    D = distance_matrix(X, quiet=quiet)
    d_c = estimate_d_c(D, fraction)
    if weights is not None:
        rho = weighed_local_density(D, d_c, weights=weights, normalize=True)
    else:
        rho = local_density(D, d_c, normalize=True)
    delta, nearest = distance_to_larger_density(D, rho, normalize=True)
    centers = cluster_centers(rho, delta, rho_min=rho_min, delta_min=delta_min, rho_x_delta_min=rho_x_delta_min)

    # Store cluster information
    for c in centers:
        clusters.append({'X': X[c,0], 'Y': X[c,1], 'rho': rho[c], 'delta': delta[c]})

    # Show scatter plot of rho vs delta
    if show_rho_vs_delta:
        plt.subplots()
        if rho_x_delta_min is not None:
            x = np.arange(0,1,0.01)
            y = np.arange(0,1,0.01)
            xx,yy = np.meshgrid(x, y)
            xy = 1.0 * ( (xx * yy) > rho_x_delta_min )
            plt.contour(x,y,xy,0, linestyles=":")
        else:
            plt.plot([rho_min,rho_min],[0,1],"r:")
            plt.plot([0,1],[delta_min,delta_min], "b:")
        plt.scatter(rho,delta)

    return clusters



def value_per_shell(X, Y, data_var, clusters, bin_size=50, start=0, end=500):
    """
    Calculates the mean value of a variable within specified distance bins from cluster centers.

    Parameters
    ----------
    X : np.array
        X-coordinates of data points.
    Y : np.array
        Y-coordinates of data points.
    data_var : np.array
        Values of the variable for each data point.
    clusters : list
        List of dictionaries with cluster information, output of the find_clusters function.
    bin_size : int, optional
        Size of the distance bins (shells), by default 50.
    start : int, optional
        Starting distance for bins, by default 0.
    end : int, optional
        Ending distance for bins, by default 500.

    Returns
    -------
    tuple
        bin_values : np.array
            Matrix with mean values of the variable per bin for each cluster [clusters, bins].
        xvalues : np.array
            Array with the center value of each bin.
    """

    # Prepare output variables
    bins = np.arange(start,end,bin_size,dtype=float)
    n_bins = len(bins)-1
    bin_values = np.zeros( (len(clusters),n_bins) )
    xvalues = bins[:-1] + (0.5*bin_size)

    # Loop over clusters
    for c_nr, c in enumerate(clusters):

        # Find distance of data points to cluster center
        D = np.sqrt( (X-c["X"])**2 + (Y-c["Y"])**2 )

        # Loop over shell-bins and calculate the mean value per shell-bin
        for b_nr in range(n_bins):
            include_ix = np.argwhere( np.logical_and( D>=bins[b_nr], D<bins[b_nr+1] ) ).ravel()
            bin_values[c_nr,b_nr] = np.nanmean(data_var[include_ix])

    # return
    return bin_values, xvalues



def distance_matrix(X, quiet=False):
    """
    Computes the distance matrix for the given data.

    Parameters
    ----------
    X : ndarray
        Data matrix X (samples x features)
    quiet : bool, optional
        Suppress print statements, by default False

    Returns
    -------
    ndarray
        Distance matrix (samples x samples)
    """

    if not quiet:
        print("Calculating pairwise distance matrix for {} samples".format(X.shape[0]))
    D = sklearn.metrics.pairwise_distances(X, metric="euclidean")
    np.fill_diagonal(D, np.NaN)
    return D



def estimate_d_c(D, fraction):
    """
    Estimates the cutoff distance d_c.

    Parameters
    ----------
    D : ndarray
        Distance matrix (samples x samples)
    fraction : float
        Fraction of data to include for distance calculation

    Returns
    -------
    float
        Cutoff distance d_c
    """

    d_array = []
    for s in range(D.shape[0]):
        d_array.append(D[s,s+1:])
    d_array = np.concatenate(d_array,axis=0)
    # print("d_array.shape={}".format(d_array.shape))
    d_c = np.percentile(d_array, fraction*100)
    # print("d_c={}".format(d_c))
    return d_c



def local_density(D, d_c, normalize=False):
    """
    Computes 'rho', the local density for each data point.

    Parameters
    ----------
    D : ndarray
        Distance matrix (samples X samples)
    d_c : float
        Cutoff distance
    normalize : bool, optional
        Normalize the local density, by default False

    Returns
    -------
    ndarray
        Local density (rho) for each data point
    """

    # Apply cuttoff
    D_cuttoff = D<d_c
    # print("D: \n{}".format(D_cuttoff[:3,:3]))

    # Some fancy weighing magic ... (see matlab script by author, they refer to this as a "Gaussian Kernel")
    rho = np.zeros((D.shape[0],))
    for s in range(len(rho)):
        rho[s] = np.sum(np.exp(-(D[s,D_cuttoff[s,:]] / d_c)**2))

    # Normalize
    if normalize:
        rho = rho / np.max(rho)

    # Calculate and return the local density vector
    return rho



def weighed_local_density(D, d_c, weights, normalize=False):
    """
    Computes a weighed version of'rho', the weighted local density for each data point.

    Parameters
    ----------
    D : ndarray
        Distance matrix (samples X samples)
    d_c : float
        Cutoff distance
    weights : ndarray
        Weights for each data point
    normalize : bool, optional
        Normalize the local density, by default False

    Returns
    -------
    ndarray
        Weighted local density (rho) for each data point
    """

    # Apply cuttoff
    D_cuttoff = D<d_c
    # print("D: \n{}".format(D_cuttoff[:3,:3]))

    # Some fancy weighing magic ... (see matlab script by author, they refer to this as a "Gaussian Kernel")
    rho = np.zeros((D.shape[0],))
    for s in range(len(rho)):

        # All distances to cells within the neighborhood d_c
        neighborhood_vector = D[s,D_cuttoff[s,:]]
        neighborhood_weights = weights[D_cuttoff[s,:]]

        # in the next line, smaller values in the neighborhood_vector result in larger values contributing to rho
        neighborhood_vector = np.exp(-(neighborhood_vector / d_c)**2)

        # Now we will scale the neighborhood by our weight
        neighborhood_vector = neighborhood_vector * neighborhood_weights

        # Finally, we sum all the scaled weights to a local weighed density
        rho[s] = np.sum(neighborhood_vector)

    # Normalize
    if normalize:
        rho = rho / np.max(rho)
        # rho = (rho-np.min(rho)) / (np.max(rho)-np.min(rho))

    # Calculate and return the local density vector
    return rho



def distance_to_larger_density(D, rho, normalize=False):
    """
    Computes 'delta', the minimum distance to a point with higher density for each data point.

    Parameters
    ----------
    D : ndarray
        Distance matrix (samples X samples)
    rho : ndarray
        Local density for each data point
    normalize : bool, optional
        Normalize the distances, by default False

    Returns
    -------
    tuple
        delta : ndarray
            Minimum distance to a point with higher density for each data point
        nearest : ndarray
            Index of the nearest point with higher density
    """

    # Output vector
    delta = np.zeros_like(rho)
    nearest = np.zeros_like(rho, dtype=np.int64)

    # Loop samples
    for s in range(len(rho)):

        # Find all samples with higher density
        nearby_ix = rho>rho[s]

        # If no samples having a larger density, set to maximum distance
        if np.sum(nearby_ix*1.0) == 0:
            delta[s] = np.nanmax(D[s,:])

        # Else set to minimum distance a point with higher density
        else:
            d_vec = D[s,:]
            d_vec[nearby_ix==False] = np.nanmax(d_vec)
            nearestby_ix = np.argmin(d_vec)
            delta[s] = D[s,nearestby_ix]
            nearest[s] = nearestby_ix

    # Normalize
    if normalize:
        delta = delta / np.max(delta)

    return delta, nearest



def cluster_centers(rho, delta, rho_min=0.2, delta_min=0.2, rho_x_delta_min=None):
    """
    Returns the indices of the cluster centers.

    Parameters
    ----------
    rho : ndarray
        Local density for each data point
    delta : ndarray
        Minimum distance to a point with higher density for each data point
    rho_min : float, optional
        Minimum rho for cluster to be included, by default 0.2
    delta_min : float, optional
        Minimum delta for cluster to be included, by default 0.2
    rho_x_delta_min : float, optional
        Minimum product of rho and delta to be included as cluster, supersedes rho_min and delta_min, by default None

    Returns
    -------
    ndarray
        Indices indicating which data points are the cluster centers
    """

    if rho_x_delta_min is None:
        centers = np.argwhere( np.logical_and( rho>rho_min, delta>delta_min ) ).ravel()
    else:
        rho_x_delta = rho*delta
        centers = np.argwhere( rho_x_delta>rho_x_delta_min ).ravel()

    center_rho_sort_ix = np.argsort(rho[centers])
    return centers[center_rho_sort_ix[::-1]]



def assign_cluster_id(rho, nearest, centers):
    """
    Assigns cluster IDs to each point.

    Parameters
    ----------
    rho : ndarray
        Local density for each data point
    nearest : ndarray
        Index of the nearest point with higher density
    centers : ndarray
        Indices of the cluster centers

    Returns
    -------
    ndarray
        Cluster IDs for each data point
    """

    order = np.argsort(rho)[::-1]
    ids = np.zeros_like(rho, dtype=np.int64)
    for s in range(len(centers)):
        ids[centers[s]] = s
    for s in range(len(rho)):
        if order[s] not in centers:
            ids[order[s]] = ids[nearest[order[s]]]
    return ids



def cluster_cores(D, d_c, rho, ids):
    """
    Returns a list of samples that are part of the cluster cores.

    Parameters
    ----------
    D : ndarray
        Distance matrix (samples X samples)
    d_c : float
        Cutoff distance
    rho : ndarray
        Local density for each data point
    ids : ndarray
        Cluster IDs for each data point

    Returns
    -------
    ndarray
        Boolean array indicating core samples
    """

    avg_border_rho = np.zeros( (len(np.unique(ids)),) )
    core = np.zeros_like(rho, dtype=bool)

    # Get matrix indicating only nearby samples by True
    D_cuttoff = D<d_c

    # Loop samples
    for s1 in range(len(rho)-1):
        for s2 in range(s1+1,len(rho)):

            # Find rho of all nearby sample without the same id (border sample)
            if ids[s1] != ids[s2] and D_cuttoff[s1,s2]:
                avg_density = 0.5 * (rho[s1] + rho[s2])
                if avg_density > avg_border_rho[ids[s1]]:
                    avg_border_rho[ids[s1]] = avg_density
                if avg_density > avg_border_rho[ids[s2]]:
                    avg_border_rho[ids[s2]] = avg_density

    for s in range(len(rho)):
        if rho[s] > avg_border_rho[ids[s]]:
            core[s] = True

    return core
