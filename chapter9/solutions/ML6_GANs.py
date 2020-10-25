# Build a Deep Convolutional Generative Adversarial Network
# Build a generator and a discriminator to generate text (in progress)

import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
# TensorFlow logger function to set verbosity & avoid warnings
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
import keras
import os
import time
from IPython import display

# Check Tensorflow version
print("Tensorflow version: {}".format(tf.__version__))

# Download datasets
print("\nDownload datasets")
mnist = keras.datasets.mnist
(train_images, train_labels), (test_images, test_labels) = mnist.load_data()
# Explore the data
print("Datasets structure:")
print("Train images: {}".format(train_images.shape))
print("Train labels: {}".format(train_labels.shape))
print("Test images: {}".format(test_images.shape))
print("Test labels: {}".format(test_labels.shape))

# Process the data before training. Scale the values to be form -1 to 1
train_images = train_images.reshape(train_images.shape[0], 28, 28, 1).astype('float32')
train_images = (train_images - 127.5) / 127.5
print("\nAfter preprocessing")
print(train_images.shape)

# Batch and suffle the data
BUFFER_SIZE = 60000
BATCH_SIZE = 256
train_dataset = tf.data.Dataset.from_tensor_slices(train_images).shuffle(BUFFER_SIZE).batch(BATCH_SIZE)
print("\nAfter batch and shuffle")
print(train_dataset)

# Build the generator model
print("\nBuilding generator")

def make_generator_model():
    
    # First layer takes a seed (random noise) as input. None is the batch size
    model = keras.Sequential()
    model.add(keras.layers.Dense(7 * 7 * 256, use_bias = False, input_shape = (100,)))
    model.add(keras.layers.BatchNormalization())
    model.add(keras.layers.LeakyReLU())

    model.add(keras.layers.Reshape((7, 7, 256)))
    assert model.output_shape == (None, 7, 7, 256) # Note: None is the batch size

    model.add(keras.layers.Conv2DTranspose(128, (5, 5), strides = (1, 1), padding = 'same', use_bias = False))
    assert model.output_shape == (None, 7, 7, 128)
    model.add(keras.layers.BatchNormalization())
    model.add(keras.layers.LeakyReLU())

    model.add(keras.layers.Conv2DTranspose(64, (5, 5), strides = (2, 2), padding = 'same', use_bias = False))
    assert model.output_shape == (None, 14, 14, 64)
    model.add(keras.layers.BatchNormalization())
    model.add(keras.layers.LeakyReLU())

    # Use as activation for the last layer the hyperbolic tangent
    model.add(keras.layers.Conv2DTranspose(1, (5, 5), strides = (2, 2), padding ='same', use_bias = False, activation = 'tanh'))
    assert model.output_shape == (None, 28, 28, 1)

    return model

# Use the (yet untrained) model to generate an image
generator = make_generator_model()
noise = tf.random.normal([1, 100])
import pdb
pdb.set_trace()
generator.trainable = False 
generated_image = generator(noise)
# Convert the Tensor to a numpy array before plotting
converted_image = keras.backend.eval(generated_image[0, :, :, 0])
plt.imshow(converted_image, cmap = "gray")
plt.show()

# Build the generator model
print("\nBuilding discriminator")

def make_discriminator_model():

    model = keras.Sequential()
    model.add(keras.layers.Conv2D(64, (5, 5), strides = (2, 2), padding = 'same', input_shape = [28, 28, 1]))
    model.add(keras.layers.LeakyReLU())
    model.add(keras.layers.Dropout(0.3))

    model.add(keras.layers.Conv2D(128, (5, 5), strides = (2, 2), padding = 'same'))
    model.add(keras.layers.LeakyReLU())
    model.add(keras.layers.Dropout(0.3))

    model.add(keras.layers.Flatten())
    model.add(keras.layers.Dense(1))

    return model

# Use the (yet untrained) discriminator to classify the generated images as real or fake
# The model will be trained to output positive values for real images, and negative values for fake images
discriminator = make_discriminator_model()
decision = discriminator(generated_image)
print("Decision: {}".format(decision))

# Discriminator Loss. Measures how well it is able to distinguish real from fakes
print("\nDiscriminator Loss")

def discriminator_loss(real_output, fake_output):

    # Compare the discriminator's predictions on real images to an array of 1s
    real_loss = keras.losses.categorical_crossentropy(tf.ones_like(real_output), real_output)
    # Compare the discriminator's predictions on fake (generated) images to an array of 0s    
    fake_loss = keras.losses.categorical_crossentropy(tf.zeros_like(fake_output), fake_output)
    total_loss = real_loss + fake_loss
    
    return total_loss

# Generator Loss. Measures how well it is able to trick the discriminator
print("\nDiscriminator Loss")

def generator_loss(fake_output):

    # Compare the decision of the discriminator with an array of ones
    return keras.losses.categorical_crossentropy(tf.ones_like(fake_output), fake_output)

# Optimizers will be different since will be train the networks separately
generator_optimizer = tf.keras.optimizers.Adam(1e-4)
discriminator_optimizer = tf.keras.optimizers.Adam(1e-4)

# Define the training loop
EPOCHS = 50
noise_dim = 100
num_examples_generate = 16
# We will reuse this seed over time to to visualize the progress in animated GIF
seed = tf.random.normal([num_examples_generate, noise_dim])

# Notice the use of tf.function, causing the function to be compiled
# @tf.function
def train_step(images):
    
    noise = tf.random.normal([BATCH_SIZE, noise_dim])

    with tf.GradientTape() as gen_tape, tf.GradientTape() as disc_tape:
        
        generator.trainable = True
        generated_images = generator(noise)

        discriminator.trainable = True
        real_output = discriminator(images)
        fake_output = discriminator(generated_images)

        gen_loss = generator_loss(fake_output)
        disc_loss = discriminator_loss(real_output, fake_output)

        # Compute the gradients of both losses
        gradients_of_generator = gen_tape.gradient(gen_loss, generator.trainable_variables)
        gradients_of_discriminator = disc_tape.gradient(disc_loss, discriminator.trainable_variables)

        # Update the parameters of both models
        generator_optimizer.apply_gradients(zip(gradients_of_generator, generator.trainable_variables))
        discriminator_optimizer.apply_gradients(zip(gradients_of_discriminator, discriminator.trainable_variables))

# Train the model
def train(dataset, epochs):

    # Iterate over the epochs recording time
    for epoch in range(epochs):
        start = time.time()

    # Iterate over batches in the dataset calling train_step()
    for image_batch in dataset:
        train_step(image_batch)

    # Produce images for the GIF as we go
    display.clear_output(wait = True)
    generate_and_save_images(generator, epoch + 1, seed)

    print ('Time for epoch {} is {} sec'.format(epoch + 1, time.time()-start))

    # Generate after the final epoch
    display.clear_output(wait = True)
    generate_and_save_images(generator, epochs, seed)

# Generate and save images
def generate_and_save_images(model, epoch, test_input):
    # Notice training is set to False.
    # This is so all layers run in inference mode (batchnorm)
    model.trainable = False
    predictions = model(test_input)
    fig = plt.figure(figsize = (4,4))

    for i in range(predictions.shape[0]):
        plt.subplot(4, 4, i+1)
        plt.imshow(predictions[i, :, :, 0] * 127.5 + 127.5, cmap='gray')
        plt.axis('off')

    plt.savefig('image_at_epoch_{:04d}.png'.format(epoch))
    plt.show()

# Train the model
print("\nTraining the model")
train(train_dataset, EPOCHS)