# Density Peak Clustering

This project implements the 'Clustering by fast search and find of density peaks' algorithm as described by Rodriguez and Laio in their paper:

> Rodriguez, Alex, and Alessandro Laio. “Clustering by Fast Search and Find of Density Peaks.” Science 344, no. 6191 (June 27, 2014): 1492–96 (https://www.science.org/doi/10.1126/science.1242072).


## Overview

Density Peak Clustering is a clustering algorithm that identifies cluster centers by finding dense regions in the data and assigns the remaining points based on their distance to these centers. This project provides a custom implementation of this algorithm. The original code supplied with the paper ("cluster_dp.m", Matlab) can be found in the demo folder. 


## Features

- Compute local density and distance to higher density points
- Identify cluster centers
- Assign cluster IDs to each point
- Determine core samples of clusters


## Demo's

The demo folder holds the following files: 
1. demo_paper_figures.ipynb: Code to reproduce a selection of figures from the original paper using this python toolbox.
2. fig1.dat: Data point coordinates of figure 1 of the paper
3. fig2_panelB.dat: Data point coordinates of figure 2B of the paper
4. fig2_panelC.dat: Data point coordinates of figure 2C of the paper
5. cluster_dp.m: Original code from the authors of the paper 

## Installation

To use this project, directly install from PyPi:

```bash
pip install densitypeakclustering
```

Or, clone the repository and install manually:

```bash
git clone https://github.com/pgoltstein/densitypeakclustering.git
cd densitypeakclustering
pip install .
```
