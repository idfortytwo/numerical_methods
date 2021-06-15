import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def exact(x):
    return (np.exp(2 - 2 * x) - 4 * np.exp(4 - 2 * x) + 4 * np.exp(2 * x) - np.exp(2 + 2 * x) - x + x * np.exp(4)) \
           / (4 - 4 * np.exp(4))


def compare_two_plots(*filenames):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.grid(linestyle='--')

    ax.set_xlabel('x')
    ax.set_ylabel('y')

    x_exact = np.linspace(0, 1, 100)
    y_exact = exact(x_exact)

    colors = ['red', 'purple']
    colors_it = iter(colors)
    markers = ['+', 'x']
    markers_it = iter(markers)
    sizes = [65, 45]
    sizes_it = iter(sizes)

    for filename in filenames:
        df = pd.read_csv(filename)
        x_header, value_header = df
        x = df[x_header]
        y_calculated = df[value_header]

        plt.scatter(x, y_calculated,
                    color=next(colors_it), marker=next(markers_it), s=next(sizes_it),
                    label=value_header)

    plt.plot(x_exact, y_exact, 'orange', label='wartość dokładna', linewidth=2)

    plt.legend(loc='best')
    # plt.show()
    plt.savefig(f'plots/{"".join(filename[:-4] for filename in filenames)}.png')


def plot_errors(filename):
    df = pd.read_csv(filename)
    print(df)

    h = df['h']
    conventional = df['conventional']
    numerov = df['Numerov']

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel('$log_{10}h$')
    ax.set_ylabel('$log_{10}|błąd|$')
    ax.grid(linestyle='--')

    plt.scatter(h, conventional, color='purple', marker='x', label='Numerov')
    plt.scatter(h, numerov, color='orange', marker='+', label='conventional')

    plt.legend(loc='best')
    # plt.show()
    plt.savefig(f'plots/{filename[:-4]}.png')


def plot_kuba(filename1, filename2):
    df1 = pd.read_csv(filename1, delimiter=' ')
    df2 = pd.read_csv(filename2, delimiter=' ')
    print(df1)
    print(df2)

    h = df1['h']
    conventional = df1['conventional']
    numerov = df2['Numerov']

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel('$log_{10}h$')
    ax.set_ylabel('$log_{10}|błąd|$')
    ax.grid(linestyle='--')

    plt.scatter(h, conventional, color='red', marker='x', label='conventional')
    plt.scatter(h, numerov, color='orange', marker='+', label='Numerov')

    plt.legend(loc='best')
    plt.show()


compare_two_plots('conventional.csv', 'numerov.csv')
# compare_plot('numerov.csv')
plot_errors('errors.csv')


# plot_kuba('BledyKonw.txt', 'BledyNumerow.txt')
# plot_kuba('c.csv', 'n.csv')
