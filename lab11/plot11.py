import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def analytic(x, t):
    return 1 + np.exp(-np.pi * np.pi * t) * np.cos(np.pi * x)


def plot(filename):
    df = pd.read_csv(filename)
    t = float(list(df)[3].split('=')[1])
    x = df['x']
    value = df['calculated']
    y_exact = df['analytic']

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel('x')
    ax.set_ylabel('U(x, t)')
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.grid(linestyle='--')

    # y_exact = analytic(x, t)

    plt.scatter(x, value, color='red', marker='.', label=filename.split('_')[0], s=10)
    plt.plot(x, y_exact, 'darkblue', label='wartość analityczna', linewidth=1.5)

    plt.legend(loc='best')
    plt.title(f'Wykresy rozwiązań numerycznych dla t = {t}')
    # plt.show()
    plt.savefig(f'plots/{filename[:-4]}.png')


def plot_max_error_t(filename):
    df = pd.read_csv(filename)

    t = df['t']
    error = df['max_error']
    # analytic = df['analytic']

    # fig = plt.figure(figsize=(16, 12))
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)

    # ax.set_xlabel('$log_{10}h$')
    # ax.set_ylabel('$log_{10}|błąd|$')
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.grid(linestyle='--')

    plt.scatter(t, error, color='red', marker='.', label=filename.split('_')[0], s=5)
    plt.legend(loc='best')
    plt.title(filename.split('_')[0])
    # plt.show()
    plt.savefig(f'plots/{filename[:-4]}.png')


def plot_all_max_errors_t(filename_KMB, filename_LU, filename_Thomas, kek):
    df_KMB = pd.read_csv(filename_KMB)
    df_LU = pd.read_csv(filename_LU)
    df_Thomas = pd.read_csv(filename_Thomas)

    n = 1
    df_KMB = df_KMB[df_KMB.index % n == 0]
    df_LU = df_LU[df_LU.index % n == 0]
    df_Thomas = df_Thomas[df_Thomas.index % n == 0]

    t_KMB = df_KMB['t']
    t_LU = df_LU['t']
    t_Thomas = df_Thomas['t']
    error_KMB = df_KMB['max_error']
    error_LU = df_LU['max_error']
    error_Thomas = df_Thomas['max_error']

    # fig = plt.figure(figsize=(16, 12))
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel('t')
    ax.set_ylabel('błąd')
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.grid(linestyle='--')

    s=10
    plt.scatter(t_KMB, error_KMB, color='slateblue', marker='+', label=f'KMB {df_KMB.columns[2]} {df_KMB.columns[3]}', s=s)
    plt.scatter(t_LU, error_LU, color='red', marker='x', label=f'LU {df_LU.columns[2]} {df_LU.columns[3]}', s=s)
    plt.scatter(t_Thomas, error_Thomas, color='gold', marker='+', label=f'Thomas {df_Thomas.columns[2]} {df_Thomas.columns[3]}', s=s)

    plt.legend(loc='best')
    plt.title("Wykres zależności maksymalnej wartości bezwzględnej błędu od czasu t")
    # plt.show()
    plt.savefig(f'plots/{kek}.png')


# def rzad(df, col, i, j):
#     calc = df.columns[col], ((df[df.columns[col]][i]) - df[df.columns[col]][j]) / (
#             df[df.columns[0]][i] - df[df.columns[0]][j])
#     print(f'{calc[0]} {calc[1]:.3f}')


def plot_max_error_h(filename):
    df = pd.read_csv(filename)

    h = df['h']
    kmb = df['kmb']
    lu = df['lu']
    thomas = df['thomas']

    # fig = plt.figure(figsize=(16, 12))
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1)

    ax.set_xlabel('$log_{10}h$')
    ax.set_ylabel('$log_{10}|błąd|$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(linestyle='--')

    plt.scatter(h, kmb, color='red', marker='x', label='KMB', s=50)
    plt.scatter(h, lu, color='slateblue', marker='+', label='LU', s=50)
    plt.scatter(h, thomas, color='orange', marker='x', label='Thomas', s=50)

    col, i, j = 1, 2, 1
    kek = ((df[df.columns[col]][i]) - df[df.columns[col]][j]) / (df[df.columns[0]][i] - df[df.columns[0]][j])
    # print(df.columns[col], kek)
    # print((np.log10(kmb[4]) - np.log10(kmb[2])) / (np.log10(h[4]) - np.log10(h[2])))

    plt.legend(loc='best')
    plt.title('Wykres zależności maksymalnej wartości bezwzględnej błędu obserwowanej dla t = t_max w funkcji kroku przestrzennego h')
    # plt.show()
    plt.savefig(f'plots/{filename[:-4]}.png')


plot_max_error_h('errors_h.csv')

# plot('Thomas_max_t.csv')
# plot_max_error_t('kmb_error_t.csv')
# plot_max_error_t('thomas_error_t.csv')
# plot_max_error_t('LU_error_t.csv')

plot_all_max_errors_t('kmb_error_t_small.csv', 'LU_error_t_small.csv', 'thomas_error_t_small.csv', 'smaller')
plot_all_max_errors_t('kmb_error_t.csv', 'LU_error_t.csv', 'thomas_error_t.csv', 'normal')

plot('kmb_1_t.csv')
plot('kmb_2_t.csv')
plot('kmb_3_t.csv')

plot('thomas_1_t.csv')
plot('thomas_2_t.csv')
plot('thomas_3_t.csv')

plot('lu_1_t.csv')
plot('lu_2_t.csv')
plot('lu_3_t.csv')