#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Script tests the clustering by density peaks algorithm on the synthetic datasets
from the paper

boundaries for fig1 set: X=[0.708, 154.567], Y=[1.133, 106.885]

Created on Sunday 1 Aug 2021

@author: pgoltstein
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import densityclustering as dc

useplot = 2

if useplot == 1:
    rho_min = 2
    delta_min = 0.2
    colors = ["#00FF00","#0000FF"]

    x = []
    y = []
    with open('fig1.dat') as f:
        for line in f:
            data = line.split()
            x.append(float(data[1]))
            y.append(float(data[2]))
    x = np.array(x)
    y = np.array(y)
    x = (x-0.708)/154.567
    y = np.abs(((y-1.133)/106.885) - 1)
    x = x * ((154.567-0.708)/(106.885-1.133))

    X = np.stack((x,y),axis=1)
    D = dc.distance_matrix(X)
    d_c = 0.2
    rho = dc.local_density(D, d_c)
    print(rho[:10])

else:
    colormap = matplotlib.cm.get_cmap('Set1')
    colors = []
    for i in range(5):
        colors.append(colormap(i/5))

    if useplot == 2:
        filename = 'fig2_panelB.dat'
        rho_min = 2
        delta_min = 0.2
    elif useplot == 3:
        filename = 'fig2_panelC.dat'
        rho_min = 20
        delta_min = 0.2
    x = []
    y = []
    with open(filename) as f:
        for line in f:
            data = line.split()
            x.append(float(data[0]))
            y.append(float(data[1]))
    x = np.array(x)
    y = np.array(y)


    X = np.stack((x,y),axis=1)
    fraction = 0.05

    D = dc.distance_matrix(X)
    d_c = dc.estimate_d_c(D,fraction)
    rho = dc.local_density(D, d_c)
    print(rho[:10])


delta, nearest = dc.distance_to_larger_density(D, rho)
print(delta[:10])
print("nearest: {}".format(nearest.dtype))
centers = dc.cluster_centers(rho, delta, rho_min=rho_min, delta_min=delta_min)
print("centers: {}".format(centers))

ids = dc.assign_cluster_id(rho, nearest, centers)
print("ids: {}".format(ids))

core = dc.core(D, d_c, rho, ids)

plt.subplots()
plt.scatter( x, y, c="#000000", s=20, alpha=1.0, edgecolors="None" )

plt.subplots()
plt.scatter( x, y, c=rho, s=20, alpha=1.0, cmap="cool", edgecolors="None" )

plt.subplots()
plt.scatter( x, y, c=delta, s=20, alpha=1.0, cmap="cool", edgecolors="None" )

plt.subplots()
# plt.scatter(x,y)
for p in range(len(ids)):
    if core[p]:
        plt.plot(x[p],y[p],marker="o",markersize=5,mfc=colors[ids[p]],mec="None")
    else:
        plt.plot(x[p],y[p],marker="o",markersize=5,mfc="#999999",mec="None")
for p in centers:
    plt.plot(x[p],y[p],marker="o",markersize=20,mfc="None",mec="#FF0000")


print(colors)

plt.subplots()
plt.scatter(rho,delta)

plt.show()
