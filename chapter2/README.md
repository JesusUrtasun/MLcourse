# Chapter 2

## Introduction to machine learning - Linear models

*Jesús Urtasun Elizari - University of Milan - 2019/20.*

Before starting we suggest to create a folder for chapter 1 where you can save all files that will be created for the exercise
```bash
cd ~/           # go the home directory
mkdir chapter2  # create the directory chapter2 in home
cd chapter2     # go inside chapter2
```


## Exercise 1 - Explore data

In this chapter we will build a basic machine learning model that performs a linear regression with some provided data.
The first exercise we will explore the input files inside the folder `data`.

1. Import `numpy` as `np` and `matplotlib.pyplot` as `plt`. 

2. Load the datasets using the `load()` method of numpy
```python
xy_training = np.load("../data/linear_training.npy")
xy_testing = np.load("../data/linear_testing.npy")
```

3. Split the training and testing sets into training and testing sets as follows.
```python
x_train = xy_training[:, 0]
y_train = xy_training[:, 1]
x_test = xy_testing[:, 0]
y_test = xy_testing[:, 1]
```

4. Plot the datasrts using `matplotlib`. 


## Exercise 2 - Prediction and loss

1. Write a loss function `MSE` that computes the minimum squared error.
```python
def MSE(y_true, y_pred):

    N = y_true.shape[0]
    delta = 
    mse = 

    return mse
```
Expected outputs: (...)

2. Write a function `predict` that computes a prediction using a linear model 
```latex 
y(x) = a * x + b
```

```python
def predict(m, b, X, Y):
   
    # Predict and compute the loss
    N = Y.shape[0]
    y_pred = X * m + b
    loss = MSE(Y, y_pred)
    
    # Compute the gradients
    dm = 
    db = 
    
    # Define a dictionary storing the gradients
    grads = {"dm": dm, "db": db}
    
    return grads, loss
```
Expected outputs: (...)


## Exercise 3: Training

1. Write a `train` function
```python
def train(m, b, X, Y, learning_rate, num_iterations = 1000, print_loss = True):

    loss_list = []
    slopes = []
    biases = []
    
    for i in range(num_iterations):
        
        # Save the initial parameter in a list before starting training
        slopes.append(m)
        biases.append(b)
        
        # Predict and compute loss
        grads, loss = 
        # Store the loss in a list after each update
        loss_list.append(loss)

        # Obtain the derivatives from the grads dictionary
        dm = grads["dm"]
        db = grads["db"]
        
        # Update rule (actual training)
        m = 
        b = 
        
        # Print the loss every 100 training iterations
        if print_loss and i % 100 == 0:
            print("Iteration = {}, Loss = {}".format(i, loss))
            
    params = {"m": m, "b": b}
    grads = {"dm": dm,"db": db}
    
    return params, grads, loss_list, slopes, biases
```
Expected outputs(...)

2. Plot the loss function with every iteration, also the evolution of the slope and bias of the linear prediction.

3. Plot the test data with the final prediction.