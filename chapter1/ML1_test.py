# ML course - Chapter 0
# Install numpy, TF and Keras

import numpy as np
import tensorflow as tf
import keras as K
# TensorFlow logger function to avoid warnings
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

# Operations with numpy
print("\n1. Operations with numpy")
a = np.random.random((2, 3))
b = np.random.random((3, 2))
c = np.dot(a, b)
print("a =\n{}\nb =\n{}\na · b =\n{}".format(a, b, c))

# Check version of tensorflow
print("\n2. Tensorflow version: {}".format(tf.__version__))

# Build a model with keras
print("\n3. Building model with Keras")
model = K.Sequential()

# End of the checks
print("\nEverything up to date")