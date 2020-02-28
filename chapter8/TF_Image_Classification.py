# Machine Learning - image classification
# Build a Neural Network using Keras (A efectos practicos, Keras y TF son lo mismo) and that performs image recognition
# Nota para acordarme: Lineas 7, 9, 12 pueden comentarse. TF solo se usa para checkear la version y mutear warnings

import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
# TensorFlow logger function to set verbosity & avoid warnings
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
import keras as K

# Check version of tensorflow
print("Tensorflow version: {}".format(tf.__version__))

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
example_input = int(input("\nChoose a picture from the training set from 0 to 60000: "))
plt.figure()
# print(train_images[example_number])
plt.imshow(train_images[example_input])
plt.colorbar()
plt.show()
# Process the data before training. Scale the values to be form 0 to 1
train_images = train_images / 255
test_images = test_images / 255

# Build the model. Specify the NN structure: number of layers, input and output shape
print("\n1. Building model. Network with one hidden layer")
model = K.Sequential([
    K.layers.Flatten(input_shape = (28, 28)),
    K.layers.Dense(128, activation = "relu"),
    K.layers.Dense(10, activation = "softmax")
    ])
model.summary()

# Compile the model. Choose the optimizer, the loss and the way training is supervised
print("\n2. Comping model")
model.compile(optimizer = "adam", loss = "sparse_categorical_crossentropy", metrics = ["accuracy"])

# Train the model. Give it the training images with their corresponding labels
print("\n3. Training model")
model.fit(train_images, train_labels, epochs = 3)

# Evaluate accuracy
print("\n4. Evaluate accuracy")
test_acc = model.evaluate(test_images, test_labels)[1]
print("Test accuracy: {}".format(test_acc))

# Making predictions. Give the test set to the model, to predict a label for every image in the set
# Predictions is a set of ten-dimensional arrays, each with numbers from 0 to 1  
predictions = model.predict(test_images)
example_test = int(input("\nChoose a picture from the testing set from 0 to 10000: "))
plt.figure()
plt.imshow(test_images[example_test])
plt.colorbar()
plt.show()
# The argmax function takes the higher value of the array
print("Prediction: {}".format(np.argmax(predictions[example_test])))
print("Truth: {}".format(test_labels[example_test]))