import pandas as pd
import matplotlib.pyplot as plt
import sys

def plotFunc(data, x_col, y_cols, xlabel, ylabel, titles):
    fig, axes = plt.subplots(2, 3, figsize=(10, 8))  # 2x2 grid of subplots
    axes = axes.flatten()  # Flatten to easily iterate

    for ax, col, title in zip(axes, y_cols, titles):
        ax.plot(data[x_col], data[col])
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    data = pd.read_csv(sys.argv[1])

    # Define column names and titles
    y_columns = ['fitness value', 'correction depth', 'counting depth', 'encoding cost', 'minimum distance', 'girth']
    titles = ['Fitness Value vs. Iteration', 'Correction Depth vs. Iteration', 'Counting Depth vs. Iteration', 
              'Encoding vs. Iteration', 'Minimum Distance vs. Iteration', 'Girth vs. Iteration', ]

    # Plot data in 2x2 subplots
    plotFunc(data, 'iteration', y_columns, 'Iteration', 'Values', titles)