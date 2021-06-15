import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# def make_plot(name):
#     df = pd.read_csv(f'{name}.csv')
#
#     fig, ax = plt.subplots(figsize=(12, 8), dpi=100)
#
#     column_colors = {
#         df.columns[1]: ('red', '+'),
#         df.columns[2]: ('orange', 'x'),
#         df.columns[3]: ('plum', '+'),
#         df.columns[4]: ('purple', 'x'),
#         df.columns[5]: ('yellowgreen', '+'),
#         df.columns[6]: ('springgreen', '+'),
#         df.columns[7]: ('gold', '+'),
#         df.columns[8]: ('dodgerblue', 'x'),
#         df.columns[9]: ('blue', '+')
#     }
#
#     for column, (color, marker) in column_colors.items():
#         df.plot(kind='scatter', x='h', y=column, color=color, marker=marker, s=25, label=column, ax=ax)
#     ax.set_yscale('log')
#     ax.set_xscale('log')
#
#     ax.yaxis.set_major_locator(ticker.LogLocator(numticks=20))
#     ax.xaxis.set_major_locator(ticker.LogLocator(numticks=20))
#
#     ax.set_title(rf'Obliczenia dla $\bf{name}$')
#     ax.set_xlabel(r'Krok siatki $h$')
#     ax.set_ylabel(r'Błąd bezwzględny przybliżeń różnicowych')
#
#     ax.grid(linestyle='--')
#
#     plt.legend(loc='best', bbox_to_anchor=(0., 0., 0.5, 0.5))
#
#     for i in range(1, 10):
#         rzad(df, i)
#
#     plt.show()


def rzad(df, col, i, j):
    calc = df.columns[col], ((df[df.columns[col]][i]) - df[df.columns[col]][j]) / (
            df[df.columns[0]][i] - df[df.columns[0]][j])
    print(f'{calc[0]} {calc[1]:.3f}')


def exact(t):
    return 1 - np.exp(-10*(t + np.arctan(t)))


def compare_plot(filename):
    df = pd.read_csv(filename)
    t_header, value_header = df
    t = df[t_header]
    y_calculated = df[value_header]

    t_exact = np.linspace(0, 3, 100)
    y_exact = exact(t_exact)

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)
    ax.grid(linestyle='--')

    ax.set_xlabel('t')
    ax.set_ylabel('y')

    plt.plot(t_exact, y_exact, 'orange', label='wartość dokładna', linewidth=2)
    plt.scatter(t, y_calculated, color='red', marker='x', s=25, label=value_header)

    plt.legend(loc='best')
    # plt.show()
    plt.savefig(f'plots/{filename[:-4]}.png')


def plot_errors(filename):
    df = pd.read_csv(filename)
    print(df)

    h = df['h']
    BME = df['BME']
    PME = df['PME']
    PMT = df['PMT']

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel('$log_{10}h$')
    ax.set_ylabel('$log_{10}|błąd|$')
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.grid(linestyle='--')

    plt.scatter(h, BME, color='red', marker='x', label='BME')
    plt.scatter(h, PME, color='orange', marker='+', label='PME')
    plt.scatter(h, PMT, color='green', marker='+', label='PMT')

    # rzad(df, 1, 0, 10)
    # rzad(df, 2, 0, 10)
    # rzad(df, 3, 0, 10)

    plt.legend(loc='best')
    # plt.show()
    plt.savefig(f'plots/{filename[:-4]}.png')


# compare_plot('PME.csv')
# compare_plot('PMT.csv')
# compare_plot('BME_stable.csv')
# compare_plot('BME_unstable.csv')
plot_errors('janek.csv')
# plot_errors('BME_PME_PMT_errors - COPY.csv')
