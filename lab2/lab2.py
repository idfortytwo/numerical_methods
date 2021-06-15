from math import exp, fabs


def f(x):
    return (1 - exp(-x)) / x


with open('cmake-build-debug/dane') as file:
    last_x = 0

    for line in file:
        x, y = line.split()[1:]
        x = float(x)
        y = float(y)
        rd = f(x)

        wzgledny = fabs(y - rd) / y

        if wzgledny > 10**(-16):
            print(f'{x}: {wzgledny}')
        else:
            last_x = x
            break

    print(f'{last_x=}')