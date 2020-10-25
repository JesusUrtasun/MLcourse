# Machine Learning - Chapter 3: Basic NN in OOP
# Build a fully operational Neural Network using OOP
# Integer vs one-hot encoding

import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow import keras as K
from sklearn import datasets
from sklearn.utils import shuffle

##### Function & Class definitions ##### 

# Activation function
def sigmoid(x):
    return 1 / (1 + np.exp(-x))

# Derivative of the activation
def der_sigmoid(x):
    return sigmoid(x) * (1 - sigmoid(x))

# Abstract class Neural Network
class Neural_Network:

    # Init method, specifying necesary inputs to instanciate the class
    def __init__(self, example_in, example_out, nlayers = 4, nodes = 128, learning = 0.25):
    
        self.input_size = len(example_in)
        self.output_size = len(example_out)
        self.nlayers = nlayers
        self.nodes = nodes
        self.learning = learning
        self.weights = []
        self.bias = []

        # Initialize the Neural Network structure. Number of layers and nodes per layer
        print("\nInitialize the Neural Network")

        # Iterate over the number of hidden layers
        hidden_layers = self.nlayers - 2
        list_nodes = [len(example_in)]
        # Hidden layers will have all the same number of nodes, given by self.nodes
        for _ in range(hidden_layers):
            input_nodes = self.nodes
            list_nodes.append(input_nodes)
        # Output layer will have the same number of nodes as the elements of the output example
        input_nodes = len(example_out)
        list_nodes.append(input_nodes)
        print("NN structure: {}".format(list_nodes))

        # Initialize the weights and biases
        for i in range(self.nlayers - 1):
            self.weights.append( np.random.rand(list_nodes[i + 1], list_nodes[i]) - 0.5 )
            self.bias.append( np.random.rand(list_nodes[i + 1]) - 0.5 )
        
    # Predict function. Also called "forward feeding". Predict from an input vector
    def predict(self, x):
        
        # Make sure x is an array
        x_array = np.array(x)
        z = [0]
        a = [x_array]
        # Iterate over the number of layers building the z, a lists with the action of each layer
        for i in range(self.nlayers - 1):
            z.append( np.dot(self.weights[i], a[-1]) + self.bias[i] )
            a.append( sigmoid(z[-1]) )
        
        return z, a

    # Loss function
    def loss(self, prediction, true):
        
        # Compute a chi2 loss function
        delta = true - prediction
        chi2 = np.dot(delta, delta) / len(prediction)
        
        return chi2

    # Backpropagation. Compute gradient of the loss and update the parameters of the network
    def backpropagation(self, z, a, true):
        
        # Gadient of the loss with respect to the output -> (true - prediction)
        delta_L = (true - a[-1])
        # Derivative of y with respect to the previous layer -> der_sigmoid(w * previous output)
        delta_y = der_sigmoid(z[-1])
        # Product delta_L and delta_y. Notice delta_y is diagonal, then the * product is valid element by element
        delta = delta_L * delta_y
        # Derivative of the layer without activation with respect to the weights w_ij -> input_vector. 
        gradient = np.outer(delta, a[-2])
        update_weights = [(delta, gradient)]
        
        # Run backwards from the output layer
        for i, weight in enumerate(reversed(self.weights[1:])):
            delta = np.dot(delta, weight) * der_sigmoid(z[-2 - i])
            gradient = np.outer(delta, a[-3 - i])
            update_weights.append((delta, gradient))
        
        # Reverse the updated list of weights and biases
        update_weights.reverse()

        # Update the lists containing the weights and biases
        for weight, bias, update in zip(self.weights, self.bias, update_weights):
            weight += self.learning * update[1]
            bias += self.learning * update[0]

        
    # Train function. Compute the loss function and backpropagates over a training set
    def train(self, input_set, target_set):
        
        # Run for each element in the input and training set
        for x, y in zip(input_set, target_set):
            z, a = self.predict(x)
            self.backpropagation(z, a, y)

    # Split the input and target sets into training and testing sets
    def train_wrapper(self, input_set, target_set, epochs, gen_test = 0.1):

        total_len = len(input_set)
        test_len = int(total_len * gen_test)

        # Test sets from the dataset. From test_length to the end
        test_x = input_set[-test_len:]
        test_y = target_set[-test_len:]

        # Train sets from the dataset. From the beginning to total - test_len
        train_x = input_set[:-test_len]
        train_y = target_set[:-test_len]

        # Shuffle the train sets before training
        for i in range(epochs):
            set_x, set_y = shuffle(train_x, train_y)
            self.train(set_x, set_y)
            print("\nTraning number {0} finished".format(i+1))
            print("Computing test")
            ave_cost, ratio = self.test(test_x, test_y)
            print("The average error was: {0}".format(ave_cost))
            print("With success ratio of: {0}".format(ratio))  

    # Test function. Define succes and compute predictions over a testing set
    def test(self, input_set, testing_set):
        
        # Define success
        loss = 0.0
        success = 0.0
        for x, y in zip(input_set, testing_set):
            _, a = self.predict(x)
            loss += self.loss(a[-1], y)
            # Accumulate one succes if prediction and example match
            if np.argmax(a[-1]) == np.argmax(y):
                success += 1.0
        
        loss = loss/len(input_set)
        succes = success/len(input_set)

        return loss, succes
    
##### Main #####

# Choose a dataset
dataset_opt = int(input("Choose a dataset\nType 1 for SKLEARN or 2 for MNIST: "))

if dataset_opt == 2:

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
    print("\nComparing raw and processed images")
    print("Train images raw: {}".format(train_images[0].shape))
    print("Train images: {}".format(train_images_input[0].shape))
    print("Test images raw: {}".format(test_images[0].shape))
    print("Test images: {}".format(test_images_input[0].shape))
    # Turn the train and test labels into one hot encoded vectors
    train_labels_hot = []
    test_labels_hot = []
    for train_label, test_label in zip(train_labels, test_labels):
        # One hot encode the train labels
        train_hot = np.zeros(10)
        train_hot[train_label] = 1
        train_labels_hot.append(train_hot)
        # One hot encode the test labels
        test_hot = np.zeros(10)
        test_hot[test_label] = 1
        test_labels_hot.append(test_hot)
    # Define a particular example to plot and check
    print("\nComparing integer and one-hot encoded labels")
    print("Train labels hot: {}".format(train_labels_hot[0]))
    print("Train labels: {}".format(train_labels[0]))
    print("Test labels hot: {}".format(test_labels_hot[0]))
    print("Test labels: {}".format(test_labels[0]))
    # Plot an example from the training set
    plt.figure()
    plt.imshow(train_images[0])
    plt.colorbar()
    plt.show()
    # For fitting in the previous architecture
    example_input = train_images_input[11]
    example_output = train_labels_hot[11]

elif dataset_opt == 1:

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
    for x, y in zip(x_set_full, y_set_full):
        x_set.append(x)
        y_set.append(y)
    # Define a particular example to plot and check
    example_number = 11
    example_input = x_set[example_number]
    example_output = y_set[example_number]
    print("\nInput {}:\n{}\nOutput {}:\n{}".format(example_input.shape, example_input, example_output.shape, example_output))
    # Reshape the input elements form 64-dim array to (8,8) matrices and plot it
    example_image = example_input.reshape(8, 8)
    plt.imshow(example_image)
    plt.colorbar()
    plt.show()

# Build the model. Instance the Neural Network class with given inputs
print("\n1. Building model")
model = Neural_Network(example_input, example_output, nlayers = 4, nodes = 128)

# Call the predict method to give a first prediction
print("\n2. First prediction")
z, a = model.predict(example_input)
# print("First z: {}\nFirst a (Input): {}".format(z[0], a[0]))
# print("Last z: {}\nLast a (Prediction): {}".format(z[-1], a[-1]))
print("Predicted output: {}".format(np.argmax(a[-1])))

# Train the model. Call the train method to run the backpropagation over the whole train set

if dataset_opt == 2:
        
    print("\n3. Training for MNIST still to be implemented")

elif dataset_opt == 1:

    print("\n3. Training model")
    # model.train(train_images, train_label)
    model.train_wrapper(x_set, y_set, epochs = 50)

    # Making predictions
    example_test =int(input("\nChoose an element to test the model: "))
    plt.figure()
    plt.imshow(x_set[example_test].reshape(8, 8))
    plt.colorbar()
    plt.show()
    z, prediction = model.predict(x_set[example_test])
    # The argmax function takes the higher value of the array
    print("\nPrediction: {}".format(np.argmax(prediction[-1])))
    print("Truth: {}".format(np.argmax(y_set[example_test])))