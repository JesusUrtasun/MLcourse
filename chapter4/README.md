# Chapter 4

### Jesús Urtasun Elizari - Uiverisity of Milan - 2019/20

First contact with TensorFlow.
Write a Python code that performs image classification with a TensorFlow model.

## Exercise 1 - Import needed libraries

Import Numpy, Matlplotlib, TensorFlow and Keras. Download and standarize the MNIST dataset.
```python
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
import keras as K

# Check version of tensorflow
print("Tensorflow version: {}".format(tf.__version__))

# Download datasets
print("\nDownload datasets")
mnist = 
(train_images, train_labels), (test_images, test_labels) = 
```
## Exercise 2

Build a Keras Sequential model with
- Input layer that converts the input (28 x 28) images into arrays of (128)
- Dense layer (128)
- Dense layer (10)

Compile it with Adam optimizer, sparse categorical cross entropy as loss, and accuracy as metrics

Train with the MNIST dataset. Check how different number of epochs affects performance.

Test accuracy

```python
# Build the model. Specify the NN structure. Number of layers, input and output shape
print("\n1. Building model. MLP with one hidden layer")
model = 

# Compile the model. Choose the optimizer, the loss and the way training is supervised
print("\n2. Comping model")
model.compile(optimizer = , loss = , metrics = )

# Train the model. Give it the training images with their corresponding labels
print("\n3. Training model")
model.fit(train_images, train_labels, epochs = )

# Evaluate accuracy
print("\n4. Evaluate accuracy")
test_acc = 
```

## Exercise 3 - Make predictions

Make predictions over the test set and ask the user to print prediction and label of one single element
```python
predictions = 
example_test = int(input("\nChoose a picture from the testing set from 0 to 10000: "))
plt.figure()
plt.imshow(test_images[example_test])
plt.colorbar()
plt.show()
```