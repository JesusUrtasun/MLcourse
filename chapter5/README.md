# Chapter 5

### Jesús Urtasun Elizari - Uiverisity of Milan - 2019/20

Build a fully operational NN using Object Oriented Programming (OOP).

## Exercise 1 - Import needed libraries

Import Numpy, Matlplotlib, TensorFlow, Keras and SKLEARN datasets.
```python
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow import keras as K
from sklearn import datasets
from sklearn.utils import shuffle
```

## Exercise 2 - OOP Neural Newtork

Build functions for the sigmoid and its derivative, then a class for the neural network.
Hint: use the following structure with a Sigmoid as activation function.

```python
# Activation function
def sigmoid(x):

# Derivative of the activation
def der_sigmoid(x):

# Abstract class Neural Network
class Neural_Network:

    # Init method, specifying necesary inputs to instanciate the class
    def __init__(self, example_in, example_out, nlayers = 4, nodes = 128, learning = 0.25):
    
        self.input_size = 
        self.output_size = 
        self.nlayers = 
        self.nodes = 
        self.learning = 
        self.weights = 
        self.bias =

    # Predict function. Also called "forward feeding". Predict from an input vector
    def predict(self, x):

    # Loss function
    def loss(self, prediction, true):

    # Backpropagation. Compute gradient of the loss and update the parameters of the network
    def backpropagation(self, z, a, true):
        
    # Train function. Compute the loss function and backpropagates over a training set
    def train(self, input_set, target_set):

    # Split the input and target sets into training and testing sets
    def train_wrapper(self, input_set, target_set, epochs, gen_test = 0.1):

    # Test function. Define succes and compute predictions over a testing set
    def test(self, input_set, testing_set):
```

## Exercise 3 - Computing gradient and backpropagation

For the backpropagation we will compute the gradient of the loss, then update the weights of the network.
```python
 # Backpropagation. Compute gradient of the loss and update the parameters of the network
    def backpropagation(self, z, a, true):
        
        # Gadient of the loss with respect to the output -> (true - prediction)
        delta_L = 
        # Derivative of y with respect to the previous layer -> der_sigmoid(w * previous output)
        delta_y = 
        # Product delta_L and delta_y. Notice delta_y is diagonal, then the * product is valid element by element
        delta = 
        # Derivative of the layer without activation with respect to the weights w_ij -> input_vector. 
        gradient = 
        update_weights = 
        
        # Run backwards from the output layer
        for i, weight in enumerate(reversed(self.weights[1:])):
            delta = 
            gradient = 
            update_weights.append((delta, gradient))
        
        # Reverse the updated list of weights and biases
        update_weights.reverse()

        # Update the lists containing the weights and biases
        for weight, bias, update in zip(self.weights, self.bias, update_weights):
            weight += self.learning * update[1]
            bias += self.learning * update[0]
```
Use the train_wrapper() function to split the dataset into a training and testing set, then shuffel before training.
```python
# Split the input and target sets into training and testing sets
    def train_wrapper(self, input_set, target_set, epochs, gen_test = 0.1):

        total_len = len(input_set)
        test_len = int(total_len * gen_test)

        # Test sets from the dataset. From test_length to the end
        test_x = 
        test_y = 

        # Train sets from the dataset. From the beginning to total - test_len
        train_x = 
        train_y = 

        # Shuffle the train sets before training
        for i in range(epochs):
            set_x, set_y = shuffle(train_x, train_y)
            self.train(set_x, set_y)
            
            ave_cost, ratio = self.test(test_x, test_y)
            print("The average error was: {0}".format(ave_cost))
            print("With success ratio of: {0}".format(ratio))  
```

## Exercise 4 - Make predictions

Make predictions over the test set and ask the user to print prediction and label of one single element
```python
print("\n3. Training model")
model.train_wrapper(x_set, y_set, epochs = 50)

# Making predictions
example_test =int(input("\nChoose an element to test the model: "))
plt.figure()
plt.imshow(x_set[example_test].reshape(8, 8))
plt.colorbar()
plt.show()
z, prediction = 
```

## Exercise 5 - Implement the MNIST dataset

Ask the user for an option (1 or 2) for running with SKLEARN dataset or the MNIST dataset, used in the previous chapter.
Careful, MNIST does not provided one-hot-encoded data, thereore the option will require longer preprocessing.
Also the training and back propagation will need to be modified.

```python
if dataset_opt == 1:

    # Download datasets
    print("\nDownload datasets")
    data = datasets.load_digits()
    x_set_full = data["data"]/16
    y_set_full = np.eye(10)[data["target"]]

    # Explore the data
    print("Datasets structure:")
    print("x set = {}\ny set = {}".format(x_set_full.shape, y_set_full.shape))
    # Generate two lists containing all x and y sets respectively
    x_set = []
    y_set = []

elif dataset_opt == 2:

    # Download datasets
    print("\nDownload datasets")
    mnist = K.datasets.mnist
    (train_images, train_labels), (test_images, test_labels) = mnist.load_data()
    
    # Explore the data
    print("Datasets structure:")
    print("Train images: {}".format(train_images.shape))
    print("Train labels: {}".format(train_labels.shape))
    print("Test images: {}".format(test_images.shape))
    print("Test labels: {}".format(test_labels.shape))
    # Process the data before training.
    # Scale the values to be form 0 to 1
    train_images = train_images / 255
    test_images = test_images / 255
    # Reshape the images to be arrays of one dimension
    train_images_input = []
    test_images_input = []
    for train_image, test_image in zip(train_images, test_images):
        train_input = train_image.reshape(784)
        train_images_input.append(train_input)
        test_input = test_image.reshape(784)
        test_images_input.append(test_input)

    # Turn the train and test labels into one hot encoded vectors
    train_labels_hot = []
    test_labels_hot = []
    for train_label, test_label in zip(train_labels, test_labels):
        # One hot encode the train labels
        train_hot = 
        train_hot[train_label] = 
        train_labels_hot.append(train_hot)
        # One hot encode the test labels
        test_hot = 
        test_hot[test_label] = 
        test_labels_hot.append(test_hot)
```