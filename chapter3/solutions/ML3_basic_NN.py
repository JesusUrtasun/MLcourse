# Machine Learning course - Chapter 1: Basics
# Build a Neural Network working with random numbers
# Load and use actual data

import numpy as np
import matplotlib.pyplot as plt
from sklearn import datasets

##### Function definitions ######

# Activation function
def sigmoid(x):
    return 1 / (1 + np.exp(-x))

# Derivative of the activation
def der_sigmoid(x):
    return sigmoid(x) * (1 - sigmoid(x))

# Neural network
def neural_network(input_vector, output_size, weights = None):
    if weights is None:
        # Generate a matrix of random weights to multiply the input by
        weights = np.random.rand(output_size, len(input_vector)) - 0.5
        prediction = np.dot(weights, input_vector)
    else:
        # Use the updated weights rather than random
        prediction = np.dot(weights, input_vector)
    # Use an activation function to make the prediction
    return sigmoid(prediction), weights

# Loss function
def loss(prediction, true):
    chi2 = 0
    for pred_i, true_i in zip(prediction, true):
        chi2 =+ (pred_i - true_i) * (pred_i - true_i)
    return chi2 / len(prediction)

# Training. Update the weights such as they minimize the loss
def train(input_set, train_set, epochs = 50):

    # Learning rate Updated_weghts = weights + learning * grad
    learning = 0.2

    # Run over the number of epochs
    for epoch in range(epochs):
        # Run over each element in the training set
        for input_vector, true_output in zip(input_set, train_set):

            # Call once the NN to have the prediction and the weights
            prediction, weights = neural_network(input_vector, len(true_output))
            
            # Compute gradient of the loss:
            # Derivative of the Loss with respect to the output y
            delta_L = (true_output - prediction)
            # Derivative of y with respect to the previous layer
            delta_y = der_sigmoid(np.dot(weights, input_vector))
            # Product delta_L and delta_y. Notice delta_y is diagonal, then the * product is valid element by element
            delta = delta_L * delta_y
            # Derivative of the layer without activation with respect to the weights w_ij -> input_vector. 
            # Outer product of delta (15) and input(10) -> Correct shape of the gradent (15, 10)
            gradient = np.outer(delta, input_vector)

            # Backpropagation. Update weights
            weights = weights + learning * gradient
            # Generate prediction with the updated weights
            new_pred = neural_network(x_set[0], y_size, weights = weights)[0]
        
        print("Epoch {}: Loss = {}".format(epoch, loss(new_pred, true_output)))

##### Main #####

# Choose a dataset
dataset_opt = int(input("Choose a dataset\nType 1 for random or 2 for SKLEARN images: "))

if dataset_opt == 1:

    # Generate Random datasets
    print("\nGenerate datasets")
    x_set = []
    y_set = []
    for i in range(10000):
        x_set.append(np.random.rand(10))
        y_set.append(np.random.rand(15))
    x_size = len(x_set[0])
    y_size = len(y_set[0])
    print("x element:\n{}:\ny element:\n{}".format(x_set[0], y_set[0]))

elif dataset_opt == 2:

    # Download datasets
    print("\nDownload datasets")
    data = datasets.load_digits()
    x_set_full = data["data"]
    y_set_full = np.eye(10)[data["target"]]
    # x_set contains 1797 arrays, of 64 elemens (input values) each
    # x_set contains 1797 arrays, of 10 elemens (output values) each
    # Explore the data
    print("Datasets structure:")
    print("x set = {}\ny set = {}".format(x_set_full.shape, y_set_full.shape))
    print("\nx element: {}".format(x_set_full[0]))
    print("y element: {}".format(y_set_full[0]))

    # Generate two lists containing all x and y sets respectively
    x_set = []
    y_set = []
    for x, y in zip(x_set_full, y_set_full):
        x_set.append(x)
        y_set.append(y)
    x_size = len(x_set[0])
    y_size = len(y_set[0])
    print("x = {}\ny = {}".format(x_set[0], y_set[0]))
    plt.figure()
    plt.imshow(x_set[0].reshape(8, 8))
    plt.colorbar()
    plt.show()

# Call the Neural Network to generate a first prediction
y_pred = neural_network(x_set[0], y_size)[0]
print("\nPrediction shape: {}\nTrue shape: {}".format(y_pred.shape, y_set[0].shape))
print("Prediction\n{}\nTrue:\n{}".format(y_pred, y_set[0]))
# Compute the loss function for this single prediction
chi2 = loss(y_pred, y_set[0])
print("\nLoss:\n{}".format(chi2))

# Call the train function with the whole input and training sets
print("\nTraining")
train(x_set, y_set)

# Test with a new prediction after the test
example_test = int(input("\nChoose an element to test the model: "))
if dataset_opt == 2:
    plt.figure()
    plt.imshow(x_set[example_test].reshape(8, 8))
    plt.colorbar()
    plt.show()
y_final = neural_network(x_set[example_test], y_size)[0]
print("\nPrediction\n{}\nTrue:\n{}".format(np.argmax(y_final), np.argmax(y_set[example_test])))