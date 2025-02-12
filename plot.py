import pandas as pd
import matplotlib.pyplot as plt
import sys

# Read the data

def plotFunc(data, col1, col2, xlabel, ylabel, title):
    plt.plot(data[col1], data[col2]) 
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()

if __name__ == '__main__':
    data = pd.read_csv(sys.argv[1])

    # Plot the data
    plotFunc(data, 'iteration', 'fitness value', 'Iteration', 'Fitness Value', 'Fitness Value vs. Iteration')
    plotFunc(data, 'iteration', 'correction depth', 'Iteration', 'Correction Depth', 'Correction Depth vs. Iteration')
    plotFunc(data, 'iteration', 'counting depth', 'Iteration', 'Counting Depth', 'Counting Depth vs. Iteration')
    plotFunc(data, 'iteration', 'minimum distance', 'Iteration', 'Minimum Distance', 'Minimum Distance vs. Iteration')