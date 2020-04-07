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

# Exercise 2: Prediction and loss

# Define the Loss, quantifying how far is the predicted value from the truth. Use the Mean Squared Error
print("\nLoss function")

def MSE(y_true, y_pred):

    N = y_true.shape[0]
    delta = y_true - y_pred
    mse = np.sum(delta * delta) / N

    return mse

# Test the loss function. Expected Output: 0.006666666666666678
y_true = np.array([1, 2, 3])
y_pred = np.array([1.1, 1.9, 3])
print("MSE = {}".format(str(MSE(y_true, y_pred))))

def predict(m, b, X, Y):
   
    # Predict and compute the loss
    N = Y.shape[0]
    y_pred = X * m + b
    loss = MSE(Y, y_pred)
    
    # Compute the gradients
    dm = 2 * np.sum(-(Y - y_pred) * X) / N
    db = 2 * np.sum(-(Y - y_pred)) / N
    
    # Define a dictionary storing the gradients
    grads = {"dm": dm, "db": db}
    
    return grads, loss

# Test the implementation. Expected Output: dm 15.7, db = 6.8, MSE = 15.444999999999999
m, b, X, Y = 1.0, 2.5, np.array([1, 2, 3, 4]), np.array([1.2, -2.3, 3.0, 4.5])
print("\nCompute forward feeding")
grads, loss = predict(m, b, X, Y)
print("dm = " + str(grads["dm"]))
print("db = " + str(grads["db"]))
print("MSE = " + str(loss))
print(loss)