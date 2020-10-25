# Chapter 3

### Jesús Urtasun Elizari - Uiverisity of Milan - 2019/20

Write the structure of a neural networks without using any libraries.
Use it to perform image classification.

## Exercise 1 - Basic structure

Build different functions for the sigmoid, its derivative, the neural network, the loss function and the training.
Remember that the derivative of the sigmoid will be needed for the computation of th gradient during training.
Hint: use the following structure with a Sigmoid as activation function.

```python
# Activation function
def sigmoid(x):

# Derivative of the activation
def der_sigmoid(x):

# Neural network
def neural_network(input_vector, output_size, weights = None):

# Loss function
def loss(prediction, true):

# Training. Update the weights such as they minimize the loss
def train(input_set, train_set, epochs = 50):
```

## Exercise 2 - Check performance for image classification

Generate a random dataset with numpy and check the performance.
Hint: Use the structure built in exercise 1 to perform classification over a random dataset.
```python

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

# Call the Neural Network to generate a first prediction
y_pred = 
print("\nPrediction shape: {}\nTrue shape: {}".format(y_pred.shape, y_set[0].shape))
print("Prediction\n{}\nTrue:\n{}".format(y_pred, y_set[0]))
# Compute the loss function for this single prediction
chi2 = 
print("\nLoss:\n{}".format(chi2))

# Call the train function with the whole input and training sets
print("\nTraining")

# Test with a new prediction after the test
example_test = int(input("\nChoose an element to test the model: "))

y_final = neural_network(x_set[example_test], y_size)[0]
print("\nPrediction\n{}\nTrue:\n{}".format(np.argmax(y_final), np.argmax(y_set[example_test])))
```

## Exercise 3 - Datasets from SKLEARN

Ask the user for an option (1 or 2) for running with random datasets or download SKLEARN handwritten digits.

```python
# Choose a dataset
dataset_opt = int(input("Choose a dataset\nType 1 for random or 2 for SKLEARN images: "))

if dataset_opt == 1:

    # Generate Random datasets
    print("\nGenerate datasets")
    (...)

elif dataset_opt == 2:

    # Download datasets
    print("\nDownload datasets")
    (...)

# Call the Neural Network to generate a first prediction
y_pred = neural_network(x_set[0], y_size)[0]
(...)
```