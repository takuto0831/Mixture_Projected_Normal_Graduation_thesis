# -*- coding: utf-8 -*-
"""
Created on Mon Dec 25 12:02:56 2017

@author: SHIO-160412-4
"""
from __future__ import print_function

from sklearn.datasets import fetch_20newsgroups
from sklearn.decomposition import TruncatedSVD
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import Normalizer
from sklearn import metrics

import numpy as np
from tabulate import tabulate

import logging
from sklearn.cluster import KMeans

import sys
sys.path.append('E:\git\spherecluster')
import spherecluster
from spherecluster import SphericalKMeans
from spherecluster import VonMisesFisherMixture

from sklearn.cluster import Means
from spherecluster import SphericalKMeans
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(seed=529)
X = np.random.uniform(low=-1, high=1, size=(300, 2))

km = KMeans(n_clusters=5, init="k-means++").fit(X)
skm = SphericalKMeans(n_clusters=5, init="k-means++").fit(X)

colors = ["r", "g", "y","m","b"]
colors = ["D", "s", "o", "v", "*"]
ax = plt.subplot(1, 2, 1, xlim=[-1.1, 1.1], ylim=[-1.1, 1.1])
for color, label in zip(colors, np.unique(km.labels_)):
   plot_data = X[km.labels_ == label]
   # plt.plot(plot_data[:, 0], plot_data[:, 1], "%s+" % color)
   plt.plot(plot_data[:, 0], plot_data[:, 1], "%s" % color,markersize=4)
#ax
#ax.set_title("k-means")

ax = plt.subplot(1, 2, 2, xlim=[-1.1, 1.1], ylim=[-1.1, 1.1])
for color, label in zip(colors, np.unique(skm.labels_)):
   plot_data = X[skm.labels_ == label]
   #plt.plot(plot_data[:, 0], plot_data[:, 1], "%s+" % color)
   plt.plot(plot_data[:, 0], plot_data[:, 1], "%s" % color,markersize=4)

#ax.set_title("skmeans")
plt.show()