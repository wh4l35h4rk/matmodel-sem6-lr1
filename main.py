import matplotlib.pyplot as plt

def exit_condition(experiment, x):
    match experiment:
        case 1:
            return True
        case 2:
            return x > 500
        case 3:
            return x > 600
        case 4:
            return x > 500


def func(t, T, T0, P, S, c, m, k, epsilon, sigma, T_max, T_min):
    Q1 = P * relay(T, T_max, T_min)
    Q2 = k * S * (T - T0)
    Q3 = epsilon * S * sigma * (T**4 - T0**4)
    return (Q1 - Q2 - Q3) / (c * m)


R = 1
def relay(T, T_max, T_min):
    if T_max == 0:
        return 1

    global R
    if T > T_max:
        R = 0
    elif T < T_min:
        R = 1
    return R


def rungeKutta(x0, y0, h, param_list, r_param_list, experiment):
    y = y0
    x_list = [x0]
    y_list = [y0]

    while 1:
        y_prev = y
        k1 = h * func(x0, y, *param_list, *r_param_list)
        k2 = h * func(x0 + 0.5 * h, y + 0.5 * k1, *param_list, *r_param_list)
        k3 = h * func(x0 + 0.5 * h, y + 0.5 * k2, *param_list, *r_param_list)
        k4 = h * func(x0 + h, y + k3, *param_list, *r_param_list)

        y = y + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)
        x0 = x0 + h

        x_list.append(x0)
        y_list.append(y)

        if (abs(y - y_prev) < 1e-5 or experiment == 4) and exit_condition(experiment, x0):
            break
    return x_list, y_list


if __name__ == '__main__':
    # МОДЕЛЬ БЕЗ ТЕРМОРЕГУЛЯТОРА
    T_relay_max = 0
    T_relay_min = 0

    # ИССЛЕДОВАНИЕ 1
    experiment = 1

    T0 = 298 # 25 градусов Цельсия
    P = 500 # 0.5 кВ
    S = 0.04
    c = 904 # алюминий
    m = 0.15
    epsilon = 1
    sigma = 5.67 * 10**(-8)
    k  = 6.7

    x0 = 0
    y = T0
    h = 0.1

    rk = rungeKutta(x0, y, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)
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
    plt.savefig('1.png', bbox_inches='tight')
    plt.show()
    plt.clf()


    # ИССЛЕДОВАНИЕ 2
    experiment = 2

    y = T0

    rk = rungeKutta(x0, y, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)
    x_list1 = rk[0][:(int)(500/h)]
    y_list1 = rk[1][:(int)(500/h)]
    T_max1 = y_list1[len(y_list1) - 1]

    P = 1000

    rk = rungeKutta(x0, y, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)
    x_list2 = rk[0]
    y_list2 = rk[1]
    T_max2 = y_list2[len(y_list2) - 1]

    S = 0.1

    rk = rungeKutta(x0, y, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)
    x_list3 = rk[0]
    y_list3 = rk[1]
    T_max3 = y_list3[len(y_list3) - 1]

    epsilon = 0.5
    P = 500
    S = 0.04

    rk = rungeKutta(x0, y, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)
    x_list4 = rk[0][:(int)(500/h)]
    y_list4 = rk[1][:(int)(500/h)]
    T_max4 = y_list4[len(y_list4) - 1]


    plt.xlabel(r'$t$, с')
    plt.ylabel(r'$T, \degree$K')

    plt.plot(x_list1, [T_max1 for i in range(len(y_list1))], '--', linewidth=1, color='grey')
    plt.plot(x_list2, [T_max2 for i in range(len(y_list2))], '--', linewidth=1, color='grey')
    plt.plot(x_list3, [T_max3 for i in range(len(y_list3))], '--', linewidth=1, color='grey')
    plt.plot(x_list4, [T_max4 for i in range(len(y_list4))], '--', linewidth=1, color='grey')

    plt.plot(x_list2, y_list2, label=r'$P = 1 кВт, S = 0.04 м^2, \epsilon = 1$')
    plt.plot(x_list4, y_list4, label=r'$P = 0.5 кВт, S = 0.04 м^2, \epsilon = 0.5$')
    plt.plot(x_list1, y_list1, label=r'$P = 0.5 кВт, S = 0.04 м^2, \epsilon = 1$')
    plt.plot(x_list3, y_list3, label=r'$P = 1 кВт, S = 0.1 м^2, \epsilon = 1$')

    plt.legend(loc='lower right')
    plt.grid()
    plt.savefig('2.png', bbox_inches='tight')
    plt.show()
    plt.clf()


    # ИССЛЕДОВАНИЕ 3
    experiment = 3

    epsilon = 1
    y = T0

    rk = rungeKutta(x0, y, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)
    x_list1 = rk[0][:(int)(600 / h)]
    y_list1 = rk[1][:(int)(600 / h)]
    T_max1 = y_list1[len(y_list1) - 1]

    c = 129 # золото

    rk = rungeKutta(x0, y, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)
    x_list2 = rk[0]
    y_list2 = rk[1]

    m = 0.25

    rk = rungeKutta(x0, y, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)
    x_list3 = rk[0]
    y_list3 = rk[1]

    c = 904

    rk = rungeKutta(x0, y, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)
    x_list4 = rk[0][:(int)(600/h)]
    y_list4 = rk[1][:(int)(600/h)]

    plt.xlabel(r'$t$, с')
    plt.ylabel(r'$T, \degree$K')

    plt.plot(x_list2, y_list2, label=r'$c = 129 \frac{Дж}{кг \cdot \degree K}\, m = 0.15 кг$')
    plt.plot(x_list3, y_list3, label=r'$c = 129 \frac{Дж}{кг \cdot \degree K}\, m = 0.25 кг$')
    plt.plot(x_list1, y_list1, label=r'$c = 904 \frac{Дж}{кг \cdot \degree K}\, m = 0.15 кг$')
    plt.plot(x_list4, y_list4, label=r'$c = 904 \frac{Дж}{кг \cdot \degree K}\, m = 0.25 кг$')

    plt.plot(x_list1, [T_max1 for i in range(len(y_list1))], '--', linewidth=1, color='grey')

    plt.legend(loc='lower right')
    plt.grid()
    plt.savefig('3.png', bbox_inches='tight')
    plt.show()
    plt.clf()


    # МОДЕЛЬ С ТЕРМОРЕГУЛЯТОРОМ
    experiment = 4

    m = 0.15

    y = T0

    rk = rungeKutta(x0, y, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)
    x_list1 = rk[0][:(int)(600 / h)]
    y_list1 = rk[1][:(int)(600 / h)]
    T_max1 = y_list1[len(y_list1) - 1]

    R = 1
    T_relay_max = 600
    T_relay_min = 500

    rk = rungeKutta(x0, y, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)
    x_list2 = rk[0]
    y_list2 = rk[1]
    T_max2 = T_relay_max
    T_min2 = T_relay_min

    R = 1
    T_relay_max = 450
    T_relay_min = 400

    rk = rungeKutta(x0, y, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)
    x_list3 = rk[0]
    y_list3 = rk[1]
    T_max3 = T_relay_max
    T_min3 = T_relay_min

    plt.xlabel(r'$t$, с')
    plt.ylabel(r'$T, \degree$K')

    plt.plot(x_list2, [T_max2 for i in range(len(y_list2))], '--', linewidth=1, color='grey')
    plt.plot(x_list2, [T_min2 for i in range(len(y_list2))], '--', linewidth=1, color='grey')

    plt.plot(x_list3, [T_max3 for i in range(len(y_list3))], '--', linewidth=1, color='grey')
    plt.plot(x_list3, [T_min3 for i in range(len(y_list3))], '--', linewidth=1, color='grey')

    plt.plot(x_list1, [T_max1 for i in range(len(y_list1))], '--', linewidth=1, color='grey')

    plt.plot(x_list2, y_list2, label=r'$T_{max} = 600 \degree K, T_{min} = 500 \degree K$')
    plt.plot(x_list3, y_list3, label=r'$T_{max} = 450 \degree K, T_{min} = 400 \degree K$')
    plt.plot(x_list1, y_list1, label='Нагреватель без терморегулятора')

    plt.legend(loc='lower right')
    plt.grid()
    plt.savefig('4.png', bbox_inches='tight')
    plt.show()
    plt.clf()