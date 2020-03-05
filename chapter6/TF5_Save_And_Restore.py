# TensorFlow Chapter 5. Save and restore models
import os
import tensorflow as tf
from tensorflow import keras

# Check the version of TensorFlow
print("TensorFlow version: {}".format(tf.__version__))
print("TensorFlow Chapter 5: Save and Restore models")

# Download the dataset and split it in training and testing set
(train_images, train_labels), (test_images, test_labels) = keras.datasets.mnist.load_data()
INPUT_SHAPE = len(train_images[0])
train_images = train_images[:1000].reshape(-1, 28 * 28) / 255
test_images = test_images[:1000].reshape(-1, 28 * 28) / 255
train_labels = train_labels[:1000]
test_labels = test_labels[:1000]
INPUT_SHAPE = len(train_images[0])

# Function declarations
def build_model():
    model = keras.Sequential([
        keras.layers.Dense(512, activation = keras.activations.relu, input_shape = (INPUT_SHAPE,)),
        keras.layers.Dropout(0.2),
        keras.layers.Dense(10, activation = keras.activations.softmax)
    ])
    model.compile(optimizer = "adam", loss = keras.losses.sparse_categorical_crossentropy, metrics = ["accuracy"])
    return model

### Main ###

# Build a model
print("\n1. Build and compile model")
model = build_model()
model.summary()
# Save checkpoints during training
print("\n2. Save checkpoints during training")
# Authomatically savep checkpints during and after the training
checkpoint_path = "training_1/cp.ckpt"
checkpoint_dir = os.path.dirname(checkpoint_path)
# Create checkpoint callback
cp_callback = tf.keras.callbacks.ModelCheckpoint(checkpoint_path, verbose = 1)
model = build_model()
# Pass callback to training
history = model.fit(train_images, train_labels, epochs = 10,
    validation_data = (test_images, test_labels), callbacks = [cp_callback])