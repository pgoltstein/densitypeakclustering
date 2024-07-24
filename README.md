# Density Peak Clustering

This project implements the 'Clustering by fast search and find of density peaks' algorithm as described by Rodriguez and Laio in their paper:

> Rodriguez, Alex, and Alessandro Laio. “Clustering by Fast Search and Find of Density Peaks.” Science 344, no. 6191 (June 27, 2014): 1492–96.


## Overview

Density Peak Clustering is a clustering algorithm that identifies cluster centers by finding dense regions in the data and assigns the remaining points based on their distance to these centers. This project provides a custom implementation of this algorithm. The original code supplied with the paper ("cluster_dp.m", Matlab) can be found in the demo folder. 


## Features

- Compute local density and distance to higher density points
- Identify cluster centers
- Assign cluster IDs to each point
- Determine core samples of clusters


## Demo's

The demo folder holds two notebooks: 
1. "demo_functionality.ipynb": A demonstration of the functionality of the python toolbox
2. "demo_paper_figures.ipynb": Code to reproduce a selection of figures from the original paper using this python toolbox.


## Installation

To use this project, clone the repository and install the required dependencies:

```bash
git clone https://github.com/pgoltstein/densitypeakclustering.git
cd densitypeakclustering
pip install .