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

# Exercise 3: Training

# Perform the optimization. Update rule for a parameter theta is theta -> theta - alpha * dtheta
print("\nTraining")

def train(m, b, X, Y, learning_rate, num_iterations = 1000, print_loss = True):

    loss_list = []
    slopes = []
    biases = []
    
    for i in range(num_iterations):
        
        # Save the initial parameter in a list before starting training
        slopes.append(m)
        biases.append(b)
        
        # Predict and compute loss
        grads, loss = predict(m, b, X, Y)
        # Store the loss in a list after each update
        loss_list.append(loss)

        # Obtain the derivatives from the grads dictionary
        dm = grads["dm"]
        db = grads["db"]
        
        # Update rule (actual training)
        m = m - learning_rate * dm
        b = b - learning_rate * db
        
        # Print the loss every 100 training iterations
        if print_loss and i % 100 == 0:
            print("Iteration = {}, Loss = {}".format(i, loss))
            
    params = {"m": m, "b": b}
    grads = {"dm": dm,"db": db}
    
    return params, grads, loss_list, slopes, biases

# Use the train function to find the best fit line, using an initial guess for m,b = 0 
# Use a learning rate of 0.01 and at least 1000 iterations of training
params, grads, loss_list, slopes, biases = train(0.1, 0.2, x_train, y_train, learning_rate = 0.01)

# Plot the results of the training
# Expected output: Should be asymptoting to values around -2 for m and 5 for b
print("\nPlotting loss function vs the iterations")
plt.figure(figsize = (9, 3))
# The left panel shows the MSE loss, which should decrease with every training step
plt.subplot(1, 3, 1)
plt.plot(loss_list)
plt.xlabel('Iteration')
plt.ylabel('Loss')
plt.yscale('log')
# Central panel shows the value of the slope m for each training step
plt.subplot(1, 3, 2)
plt.plot(slopes)
plt.xlabel('Iteration')
plt.ylabel('Slope')
# Central panel shows the value of the bias b for each training step
plt.subplot(1, 3, 3)
plt.plot(biases)
plt.xlabel('Iteration')
plt.ylabel('Bias')
plt.tight_layout()
plt.show()

# Plot the final fit for the dta
print("\nPlotting the best fit for the data")
plt.figure()
xrange = np.linspace(0, 10)
y_pred = xrange * params['m'] + params['b']
plt.plot(xrange, y_pred)
plt.scatter(x_test, y_test, color = 'C1', label = 'testing data')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.legend()
plt.show()