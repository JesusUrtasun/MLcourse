# Chapter 1: Linear models
# Author: Jesús Urtasun - 2020

# Exercise 1: Load data

import numpy as np
import matplotlib.pyplot as plt

# Load datasets
print("\nLoad datasets")
xy_training = np.load("../data/linear_training.npy")
xy_testing = np.load("../data/linear_testing.npy")
print("train set:\n{}\ntest set:\n{}".format(xy_training, xy_testing))

# Split the dataset into a training and a testing set
x_train = xy_training[:, 0]
y_train = xy_training[:, 1]
x_test = xy_testing[:, 0]
y_test = xy_testing[:, 1]
print("train set:\n{}\ntrain labels:\n{}".format(x_train, y_train))
print("test set:\n{}\ntest labels:\n{}".format(x_test, y_test))

# Plot the data to fit once the training is done (test set). 
plt.figure()
plt.scatter(x_train, y_train, label = "training data")
plt.scatter(x_test, y_test, label = "testing data")
plt.xlabel("$x$")
plt.ylabel("$y(x)$")
plt.legend()
plt.show()