# Build a Deep Convolutional Generative Adversarial Network using keras
# Generate images to fake the MNIST dataset

import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
import keras
from keras import backend as K

# Check Tensorflow version
print("Tensorflow version: {}".format(tf.__version__))

# Function to sample from
def f(x):
    return np.dot(x, x)

# Plot a given function
def plot(function):
    y = [f(x) for x in x]
    plt.plot(x, y)
    plt.show()

# Generate a random sampling and compute their images
def generate_samples(n = 100):
    X1 = np.random.rand(n) - 0.5
    Y1 = X1 * X1
    # Stack arrays
    X1 = X1.reshape(n, 1)
    Y1 = Y1.reshape(n, 1)
    return np.hstack((X1, Y1))

# Generate real samples labeled with 1. Outputs the points and an array on ones
def generate_real_samples(n):
    # Generate inputs in [-0.5, 0.5]
    X1 = np.random.rand(n) - 0.5
    Y1 = X1 * X1
    # Stack arrays
    X1 = X1.reshape(n, 1)
    Y1 = Y1.reshape(n, 1)
    X = np.hstack((X1, Y1))
    # Generate class labels
    y = np.ones((n, 1))
    return X, y

# # Generate fake samples labeled with 0
# def generate_fake_samples(generator, n):
#     # Generate inputs in [-1, 1]
#     X1 = -1 + np.random.rand(n) * 2
#     # Generate outputs in [-1, 1]
#     X2 = -1 + np.random.rand(n) * 2
#     # Stack arrays
#     X1 = X1.reshape(n, 1)
#     X2 = X2.reshape(n, 1)
#     X = np.hstack((X1, X2))
#     # Generate class labels
#     y = np.zeros((n, 1))
#     return X, y

# Generate fake samples using the latent space
def generate_fake_samples(generator, latent_dim, n):
    
    # Generate points in latent space
    x_input = generate_latent_points(latent_dim, n)
    # Predicted outputs
    X = generator.predict(x_input)
    # Plot the results
    plt.scatter(X[:, 0], X[:, 1])
    plt.show()

# Generate points in the latent space. Array of random numbers from a standard Gaussian
def generate_latent_points(latent_dim, n):

    # Generate points in the latent space
    x_input = np.random.randn(latent_dim * n)
    # Reshape into a batch of inputs for the network
    x_input = x_input.reshape(n, latent_dim)

    return x_input

# Build generator. Take as input a point from the domain and outputs a new sample, meaning a vector of x and its image
# Define a latent space of five dimensions (...) Generate inputs by random sampling from a Gaussian distribution. Latent space of 5 dim
def build_generator(seed, n_outputs = 2):

    # Hidden layer
    model = keras.Sequential()
    model.add(keras.layers.Dense(15, activation = "relu", kernel_initializer = "he_uniform", input_dim = seed))
    model.add(keras.layers.Dense(n_outputs, activation = "linear"))
    
    return model

# Build discriminator. Take a vector of two elements and output a classification prediction real/fake
def build_discriminator(n_inputs = 2):

    model = keras.Sequential()
    model.add(keras.layers.Dense(25, activation = "relu", kernel_initializer = "he_uniform", input_dim = n_inputs))
    # model.add(keras.layers.LeakyReLU())
    model.add(keras.layers.Dense(1, activation = "sigmoid"))
    # model.add(keras.layers.LeakyReLU())
    
    # Compile model
    model.compile(loss = "binary_crossentropy", optimizer = "adam", metrics = ["accuracy"])

    return model

# Train discriminator
def train_discriminator(model, epochs = 1000, batch_size = 512):

    half_batch = int(batch_size * 0.5)
    # Run epochs manually
    for i in range(epochs):
        # Generate real examples
        X_real, y_real = generate_real_samples(half_batch)
        # Update model
        model.train_on_batch(X_real, y_real)
        # Generate fake examples
        X_fake, y_fake = generate_real_samples(half_batch)
        # Update model
        model.train_on_batch(X_fake, y_fake)
        # Evaluate model
        loss_real, acc_real = model.evaluate(X_real, y_real, verbose = 0)
        loss_fake, acc_fake = model.evaluate(X_fake, y_fake, verbose = 0)
        print("Epoch: {}, Real accuracy: {}, Fake accuracy: {}".format(i, acc_real, acc_fake))

##### Main #####

# Plot the function to learn
num_points = 10
x = []
for i in range(num_points):
    x.append(i)
plot(f(x))

# Generate samples
data = generate_samples()
# Plot samples
plt.scatter(data[:, 0], data[:, 1])
plt.show()

# Latent dimension
latent_dim = 5
# Build generator and discriminator models
print("\nBuilding generator")
generator = build_generator(latent_dim)
generator.summary()

# Generate and plot generated samples
generate_fake_samples(generator, latent_dim, 100)

print("\nBuilding discriminator")
discriminator = build_discriminator()
discriminator.summary()

# Train the discriminator
print("\nTraining discriminator")
train_discriminator(discriminator)