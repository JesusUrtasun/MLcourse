# Machine Learning - image classification
# Build a Neural Network using Keras and that performs binary classification
# Nota para acordarme: Lineas 9, 11, 15 pueden comentarse. TF solo se usa para checkear la version y mutear warnings

import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import load_breast_cancer
from sklearn.preprocessing import StandardScaler
import tensorflow as tf
# TensorFlow logger function to set verbosity & avoid warnings
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
import keras as K

# Check version of tensorflow
print("Tensorflow version: {}".format(tf.__version__))

# Download datasets
print("\nDownload datasets")
# Data is loaded as a sklearn bunch, dictionary-type object
data = load_breast_cancer()
sc = StandardScaler()
data_files = data["data"] # ya numpy array =D
data_labels = data["target"]

# Split the files and labels in to training and testing lists
# Normalize the files using the StandardScaler from sklearn
print("\nPreprocessing datasets")
train_files = sc.fit_transform(data_files[0:500])
train_labels = data_labels[0:500]
test_files = sc.fit_transform(data_files[500:569])
test_labels = data_labels[500:569]
print("Train files: {}".format(len(train_files)))
print("Train labels: {}".format(len(train_labels)))
print("Test files: {}".format(len(test_files)))
print("Test labels: {}".format(len(test_labels)))
example_input = int(input("\nChoose an element from the training set from 0 to 900: "))
print("train_files - example: {}".format(train_files[example_input]))

# Build the model. Specify the NN structure: number of layers, input and output shape
print("\n1. Building model. Network with one hidden layer")
model = K.Sequential([
    K.layers.Dense(32, activation = "relu"),
    K.layers.Dense(8, activation = "relu"),
    K.layers.Dense(1, activation = "sigmoid")
    ])
model.predict(train_files[0].reshape(1, 30))
model.summary()

# Compile the model. Choose the optimizer, the loss and the way training is supervised
print("\n2. Comping model")
model.compile(optimizer = "adam", loss = "binary_crossentropy", metrics = ["accuracy"])
print("\n2.1 Evaluate accuracy before training")
test_acc = model.evaluate(test_files, test_labels)[1]
print("Accuracy before training: {}".format(test_acc))

# Train the model. Give it the training images with their corresponding labels
print("\n3. Training model")
model.fit(train_files, train_labels, epochs = 30)

# Evaluate accuracy
print("\n4. Evaluate accuracy")
test_acc = model.evaluate(test_files, test_labels)[1]
print("Test accuracy: {}".format(test_acc))

# Making predictions. Give the test set to the model, to predict a label for every image in the set
# Predictions is a set of ten-dimensional arrays, each with numbers from 0 to 1  
predictions = model.predict(test_files)
example_test = int(input("\nChoose a file from the test set from 0 to 69: "))

if predictions[example_test] > 0.5:
    print("Prediction: {} -  Positive".format(predictions[example_test]))
else:
    print("Prediction: {} -  Negative".format(predictions[example_test]))
if test_labels[example_test] > 0.5:
    print("Truth: {} - Positive".format(test_labels[example_test]))
else:
    print("Truth: {} - Negative".format(test_labels[example_test]))