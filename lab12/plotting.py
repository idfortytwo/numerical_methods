import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def exact_value(x):
    return 1 / (1 + 10 * x ** 6)


def compare_plot(n):
    df = pd.read_csv(f'n{n}.csv')
    df_nodes_equi = pd.read_csv(f'n{n}_nodes_equi.csv')
    df_nodes_chebyshev = pd.read_csv(f'n{n}_nodes_chebyshev.csv')

    x_h, equi_h, chebyshev_h = df
    x, equi, chebyshev = df[x_h], df[equi_h], df[chebyshev_h]

    x_nodes_equi, y_nodes_equi = df_nodes_equi['x'], df_nodes_equi['y']
    x_nodes_chebyshev, y_nodes_chebyshev = df_nodes_chebyshev['x'], df_nodes_chebyshev['y']

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(f'{n} nodes')
    ax.set_xlim(-1.1, 1.1)

    x_exact = np.linspace(-1, 1, 100)
    y_exact = exact_value(x_exact)
    plt.plot(x_exact, y_exact, 'darkviolet', label='f', linewidth=1.5)

    plt.scatter(x, equi, color='orange', marker='.', s=25, label=equi_h)
    plt.scatter(x_nodes_equi, y_nodes_equi, color='orange', marker='x', s=100, label='equidistant nodes')

    plt.scatter(x, chebyshev, color='red', marker='.', s=25, label=chebyshev_h)
    plt.scatter(x_nodes_chebyshev, y_nodes_chebyshev, color='red', marker='x', s=100, label='Chebyshev nodes')

    plt.legend(loc='best')
    plt.show()
    # plt.savefig(f'plots/n{n}.png')


# compare_plot(5)
# compare_plot(10)
# compare_plot(12)
compare_plot(14)
# compare_plot(25)
