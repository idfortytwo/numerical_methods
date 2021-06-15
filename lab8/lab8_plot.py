import matplotlib.pyplot as plt
from matplotlib import ticker
import pandas as pd
import math


def rzad(df, i):
    calc = df.columns[i], ((math.log10(df[df.columns[i]][0])) - math.log10(df[df.columns[i]][5])) / (
            math.log10(df[df.columns[0]][0]) - math.log10(df[df.columns[0]][5]))
    print(f'{calc[0]} {calc[1]:.3f}')


def make_plot(name):
    df = pd.read_csv(f'{name}.csv')

    fig, ax = plt.subplots(figsize=(12, 8), dpi=100)

    column_colors = {
        df.columns[1]: ('red', '+'),
        df.columns[2]: ('orange', 'x'),
        df.columns[3]: ('plum', '+'),
        df.columns[4]: ('purple', 'x'),
        df.columns[5]: ('yellowgreen', '+'),
        df.columns[6]: ('springgreen', '+'),
        df.columns[7]: ('gold', '+'),
        df.columns[8]: ('dodgerblue', 'x'),
        df.columns[9]: ('blue', '+')
    }

    for column, (color, marker) in column_colors.items():
        df.plot(kind='scatter', x='h', y=column, color=color, marker=marker, s=25, label=column, ax=ax)
    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.yaxis.set_major_locator(ticker.LogLocator(numticks=20))
    ax.xaxis.set_major_locator(ticker.LogLocator(numticks=20))

    ax.set_title(rf'Obliczenia dla $\bf{name}$')
    ax.set_xlabel(r'Krok siatki $h$')
    ax.set_ylabel(r'Błąd bezwzględny przybliżeń różnicowych')

    ax.grid(linestyle='--')

    plt.legend(loc='best', bbox_to_anchor=(0., 0., 0.5, 0.5))

    for i in range(1, 10):
        rzad(df, i)

    plt.show()


print('double')
make_plot('double')

print('\nfloat')
make_plot('float')


