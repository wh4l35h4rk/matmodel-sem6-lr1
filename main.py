import matplotlib.pyplot as plt

def func(t, T, T0, P, S, c, m, k, epsilon, sigma):
    Q1 = P
    Q2 = k * S * (T - T0)
    Q3 = epsilon * S * sigma * (T**4 - T0**4)
    return (Q1 - Q2 - Q3) / (c * m)


def rungeKutta(x0, y0, h, param_list):
    y = y0
    x_list = [x0]
    y_list = [y0]

    while 1:
        y_prev = y
        k1 = h * func(x0, y, *param_list)
        k2 = h * func(x0 + 0.5 * h, y + 0.5 * k1, *param_list)
        k3 = h * func(x0 + 0.5 * h, y + 0.5 * k2, *param_list)
        k4 = h * func(x0 + h, y + k3, *param_list)

        y = y + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
        x0 = x0 + h

        x_list.append(x0)
        y_list.append(y)

        # if abs(y - y_prev) < 1e-5:
        if abs(y - y_prev) < 1e-5 and x0 > 500:
                break
    return x_list, y_list


def find_max(T0, P, S, k, epsilon, sigma):
    # (epsilon * S * sigma) * T ^ 4 + (k * S) * T - (k * S * T0 + epsilon * S * sigma * T0 ^ 4 + P) = 0

    a = epsilon * S * sigma
    b = 0
    c = 0
    d = k * S
    e = k * S * T0 + epsilon * S * sigma * T0**4 + P

    # a * x^4 + b * x^3 + c * x^2 + d * x = e
    # x^4 + b/a * x^3 + c/a * x^2 + d/a * x - e = 0

    p = c / a
    q = d / a
    r = - e

    # y^3 - 2p * y^2 + (p^2 - 4r) * y + q^2 = 0

    a = 1
    b = 2 * p
    c = p**2 - 4 * r
    d = q**2

    p = (3 * a * c - b**2) / (3 * a**2)
    q = (2 * b**3 - 9 * a * b * c + 27 * a**2 * d) / (27 * a**3)

    Q = (p / 3)**3 + (q / 2)**2

    alpha = (- q / 2 + Q**0.5)**(1/3)
    beta = (- q / 2 - Q**0.5)**(1/3)

    y1 = alpha + beta
    y2 = complex(-(alpha + beta) / 2,  (alpha + beta) / 2 * 3**0.5)
    y3 = complex(-(alpha + beta) / 2,  -(alpha + beta) / 2 * 3**0.5)
    return y1, y2, y3


if __name__ == '__main__':
    # ИССЛЕДОВАНИЕ 1

    T0 = 298 # 25 градусов цельсия
    P = 500 # 0.5 кВ
    S = 0.04
    c = 904 # алюминий
    m = 0.15
    epsilon = 1
    sigma = 5.67 * 10**(-8)
    k = 6.706

    x0 = 0
    y = T0
    h = 0.1

    rk = rungeKutta(x0, y, h, [T0, P, S, c, m, k, epsilon, sigma])
    x_list = rk[0]
    y_list = rk[1]
    T_max = y_list[len(y_list) - 1]
    # print(T_max - 273)

    plt.plot(x_list, y_list, label=r'$T(t, T_0, P, S, c, m, k, \epsilon, \sigma)$')
    plt.plot(x_list, [T_max for i in range(len(y_list))], '--', linewidth=1, label=r'$T_{max}$')
    plt.xlabel(r'$t$, с')
    plt.ylabel(r'$T, \degree$K')
    old_ticks = list(plt.yticks()[0])
    old_ticks.remove(650)
    plt.yticks(old_ticks + [T_max])
    plt.legend(loc='lower right')
    plt.grid()
    plt.savefig('1.png')
    plt.show()
    plt.clf()


    # ИССЛЕДОВАНИЕ 2

    T0 = 298 # 25 градусов цельсия
    P = 500 # 0.5 кВ
    S = 0.04
    c = 904 # алюминий
    m = 0.15
    epsilon = 1
    sigma = 5.67 * 10**(-8)
    k = 6.706

    y = T0

    rk = rungeKutta(x0, y, h, [T0, P, S, c, m, k, epsilon, sigma])
    x_list1 = rk[0][:(int)(500/h)]
    y_list1 = rk[1][:(int)(500/h)]
    T_max1 = y_list1[len(y_list1) - 1]

    P = 1000

    rk = rungeKutta(x0, y, h, [T0, P, S, c, m, k, epsilon, sigma])
    x_list2 = rk[0]
    y_list2 = rk[1]
    T_max2 = y_list2[len(y_list2) - 1]

    S = 0.1

    rk = rungeKutta(x0, y, h, [T0, P, S, c, m, k, epsilon, sigma])
    x_list3 = rk[0]
    y_list3 = rk[1]
    T_max3 = y_list3[len(y_list3) - 1]

    epsilon = 0.5
    P = 500
    S = 0.04

    rk = rungeKutta(x0, y, h, [T0, P, S, c, m, k, epsilon, sigma])
    x_list4 = rk[0][:(int)(500/h)]
    y_list4 = rk[1][:(int)(500/h)]
    T_max4 = y_list4[len(y_list4) - 1]


    plt.xlabel(r'$t$, с')
    plt.ylabel(r'$T, \degree$K')

    plt.plot(x_list1, y_list1, label=r'$P = 0.5 кВт, S = 0.04 м^2, \epsilon = 1$')
    plt.plot(x_list2, y_list2, label=r'$P = 1 кВт, S = 0.04 м^2, \epsilon = 1$')
    plt.plot(x_list3, y_list3, label=r'$P = 1 кВт, S = 0.1 м^2, \epsilon = 1$')
    plt.plot(x_list4, y_list4, label=r'$P = 0.5 кВт, S = 0.04 м^2, \epsilon = 0.5$')

    plt.plot(x_list1, [T_max1 for i in range(len(y_list1))], '--', linewidth=1)
    plt.plot(x_list2, [T_max2 for i in range(len(y_list2))], '--', linewidth=1)
    plt.plot(x_list3, [T_max3 for i in range(len(y_list3))], '--', linewidth=1)
    plt.plot(x_list4, [T_max4 for i in range(len(y_list4))], '--', linewidth=1)

    plt.legend(loc='lower right')
    plt.grid()
    plt.savefig('2.png')
    plt.show()


