import matplotlib.pyplot as plt

def exit_condition(experiment):
    match experiment:
        case 1:
            return 0
        case 2:
            return 500
        case 3:
            return 600
        case 4:
            return 500


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

        if (abs(y - y_prev) < 1e-5 or experiment == 4) and x0 > exit_condition(experiment):
            break
    return x_list, y_list


def exp_data(x0, h, param_list, r_param_list, experiment):
    y = T0
    rk = rungeKutta(x0, y, h, param_list, r_param_list, experiment)

    if experiment == 1:
        x_list = rk[0]
        y_list = rk[1]
        T_max = y_list[len(y_list) - 1]
        T_min = 0
    else:
        limit = exit_condition(experiment)
        x_list = rk[0][:(int)(limit / h)]
        y_list = rk[1][:(int)(limit / h)]
        if experiment == 4:
            T_max = T_relay_max
            T_min = T_relay_min
        else:
            T_max = y_list[len(y_list) - 1]
            T_min = 0

    return {"x": x_list,
            "y": y_list,
            "T_max": T_max,
            "T_min": T_min}


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

    data = exp_data(x0, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)

    plt.plot(data["x"], data["y"], label=r'$T(t, T_0, P, S, c, m, k, \epsilon, \sigma)$')
    plt.plot(data["x"], [data["T_max"] for i in range(len(data["y"]))], '--', linewidth=1, label=r'$T_{max}$')

    plt.xlabel(r'$t$, с')
    plt.ylabel(r'$T, \degree$K')
    old_ticks = list(plt.yticks()[0])
    old_ticks.remove(650)
    plt.yticks(old_ticks + [data["T_max"]])
    plt.legend(loc='lower right')
    plt.grid()
    plt.savefig('1.png', bbox_inches='tight')
    plt.show()
    plt.clf()


    # ИССЛЕДОВАНИЕ 2
    experiment = 2

    data_1 = exp_data(x0, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)

    P = 1000
    data_2 = exp_data(x0, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)

    S = 0.1
    data_3 = exp_data(x0, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)

    epsilon = 0.5
    P = 500
    S = 0.04
    data_4 = exp_data(x0, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)


    plt.xlabel(r'$t$, с')
    plt.ylabel(r'$T, \degree$K')

    plt.plot(data_1["x"], [data_1["T_max"] for i in range(len(data_1["y"]))], '--', linewidth=1, color='grey')
    plt.plot(data_2["x"], [data_2["T_max"] for i in range(len(data_2["y"]))], '--', linewidth=1, color='grey')
    plt.plot(data_3["x"], [data_3["T_max"] for i in range(len(data_3["y"]))], '--', linewidth=1, color='grey')
    plt.plot(data_4["x"], [data_4["T_max"] for i in range(len(data_4["y"]))], '--', linewidth=1, color='grey')

    plt.plot(data_2["x"], data_2["y"], label=r'$P = 1 кВт, S = 0.04 м^2, \epsilon = 1$')
    plt.plot(data_4["x"], data_4["y"], label=r'$P = 0.5 кВт, S = 0.04 м^2, \epsilon = 0.5$')
    plt.plot(data_1["x"], data_1["y"], label=r'$P = 0.5 кВт, S = 0.04 м^2, \epsilon = 1$')
    plt.plot(data_3["x"], data_3["y"], label=r'$P = 1 кВт, S = 0.1 м^2, \epsilon = 1$')

    plt.legend(loc='lower right')
    plt.grid()
    plt.savefig('2.png', bbox_inches='tight')
    plt.show()
    plt.clf()


    # ИССЛЕДОВАНИЕ 3
    experiment = 3

    epsilon = 1
    data_1 = exp_data(x0, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)

    c = 129 # золото
    data_2 = exp_data(x0, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)

    m = 0.25
    data_3 = exp_data(x0, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)

    c = 904
    data_4 = exp_data(x0, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)

    plt.xlabel(r'$t$, с')
    plt.ylabel(r'$T, \degree$K')

    plt.plot(data_2["x"], data_2["y"], label=r'$c = 129 \frac{Дж}{кг \cdot \degree K}\, m = 0.15 кг$')
    plt.plot(data_3["x"], data_3["y"], label=r'$c = 129 \frac{Дж}{кг \cdot \degree K}\, m = 0.25 кг$')
    plt.plot(data_1["x"], data_1["y"], label=r'$c = 904 \frac{Дж}{кг \cdot \degree K}\, m = 0.15 кг$')
    plt.plot(data_4["x"], data_4["y"], label=r'$c = 904 \frac{Дж}{кг \cdot \degree K}\, m = 0.25 кг$')

    plt.plot(data_1["x"], [data_1["T_max"] for i in range(len(data_1["y"]))], '--', linewidth=1, color='grey')

    plt.legend(loc='lower right')
    plt.grid()
    plt.savefig('3.png', bbox_inches='tight')
    plt.show()
    plt.clf()


    # МОДЕЛЬ С ТЕРМОРЕГУЛЯТОРОМ
    experiment = 4

    m = 0.15
    data_1 = exp_data(x0, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], 2)

    R = 1
    T_relay_max = 600
    T_relay_min = 500
    data_2 = exp_data(x0, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)

    R = 1
    T_relay_max = 450
    T_relay_min = 400
    data_3 = exp_data(x0, h, [T0, P, S, c, m, k, epsilon, sigma], [T_relay_max, T_relay_min], experiment)

    plt.xlabel(r'$t$, с')
    plt.ylabel(r'$T, \degree$K')

    plt.plot(data_2["x"], [data_2["T_max"] for i in range(len(data_2["y"]))], '--', linewidth=1, color='grey')
    plt.plot(data_2["x"], [data_2["T_min"] for i in range(len(data_2["y"]))], '--', linewidth=1, color='grey')

    plt.plot(data_3["x"], [data_3["T_max"] for i in range(len(data_3["y"]))], '--', linewidth=1, color='grey')
    plt.plot(data_3["x"], [data_3["T_min"] for i in range(len(data_3["y"]))], '--', linewidth=1, color='grey')

    plt.plot(data_1["x"], [data_1["T_max"] for i in range(len(data_1["y"]))], '--', linewidth=1, color='grey')

    plt.plot(data_2["x"], data_2["y"], label=r'$T_{max} = 600 \degree K, T_{min} = 500 \degree K$')
    plt.plot(data_3["x"], data_3["y"], label=r'$T_{max} = 450 \degree K, T_{min} = 400 \degree K$')
    plt.plot(data_1["x"], data_1["y"], label='Нагреватель без терморегулятора')

    plt.legend(loc='lower right')
    plt.grid()
    plt.savefig('4.png', bbox_inches='tight')
    plt.show()
    plt.clf()