import matplotlib.pyplot as plt
import matplotlib.patches as ptc
import matplotlib.path as pth
import numpy as np
from scipy.integrate import quad


fig, ax = plt.subplots(2,3, figsize = (30,12))

ax = [item for sublist in ax for item in sublist]


print("Сколько наборов данных?\n>")
nsets = int(input())

print("Вводите данные через пробел:\nКритическая температура, температура базы, сопротивление контакта, сопротивление сверхпроводящих линий, гамма.\nНаборы разделять Enter\n")
vars = [(float(x) for x in input().split()) for n in range(nsets)] #


def DrawGraphs(T_C, T_b, R, R_lead, gamma, num):

    k_B = 1.380649 * 10 ** (-23)
    e_0 = 1.602176 * 10 ** (-19)
    R_g = R - R_lead
    delta_0 = 1.764 * k_B * T_C

    # Сверхпроводящая щель

    def delta(temp):
        if temp >= T_C:
            return 0
        else:
            return 1.76 * k_B * T_C * np.tanh(1.74 * np.sqrt(T_C / temp - 1))

    temps = np.arange(0.0,0.7,0.002)
    ax[0].plot(temps, list(map(lambda x: delta(x) / delta_0, temps)), label = "set "+str(num))
    ax[0].set_title("Сверхпроводящая щель")
    ax[0].set_xlabel("Temp")
    ax[0].set_ylabel(r'$\frac{\Delta(Temp)}{\Delta_0}$')
    ax[0].legend(loc = 4)
    # Зависимость температуры от напряжения

    def integr(x):
        return x**2 / np.cosh(x) ** 2

    def Int(temp):
        return quad(integr, delta(temp) / 2 * k_B * temp, 200)[0]

    def subInt(temp):
        return temp * Int(temp)

    def Volt(temp):
        return 1 / 2 * (- gamma * delta(temp) / e_0 + np.sqrt((gamma * delta(temp) / e_0) ** 2 + 32 * R / R_lead * k_B ** 2 / e_0 ** 2 * quad(subInt, T_b, temp)[0]))

    temps2 = np.arange(0.05, 0.7, 0.002)
    ax[1].plot(list(map(lambda x: e_0 * Volt(x) / delta_0,temps2)), temps2, label = "set "+str(num))
    ax[1].set_title("Зависимость температуры от напряжения")
    ax[1].set_xlabel(r'$\frac{e_0 * Volt(temp)}{\Delta_0}$')
    ax[1].set_ylabel("Temp")
    ax[1].legend(loc = 4)

    ax[2].plot(list(map(lambda x: Volt(x),temps2)), temps2, label = "set "+str(num))
    ax[2].set_title("Зависимость температуры от напряжения")
    ax[2].set_xlabel("Volt(Temp)")
    ax[2].set_ylabel("Temp")
    ax[2].legend(loc = 4)
    # Зависимость тока от напряжения

    def Current(temp):
        return Volt(temp) / R + gamma * delta(temp) / (e_0 * R)

    def Ohms(temp):
        return Volt(temp) / R

    temps2 = np.arange(0.05, 0.7, 0.002)
    ax[3].plot(list(map(lambda x: Volt(x),temps2)), list(map(lambda x: Current(x),temps2)), label = "Current(Temp) "+str(num))
    ax[3].plot(list(map(lambda x: Volt(x),temps2)), list(map(lambda x: Ohms(x),temps2)), label = "Ohms(Temp) "+str(num))

    ax[3].set_title("Зависимость тока от напряжения")
    ax[3].set_xlabel("Volt(Temp)")
    ax[3].set_ylabel("I(V)")

    ax[3].ticklabel_format(style = 'scientific', scilimits = [-1, 1])

    ax[3].legend(loc = 4)
    # Проводимость от температуры


    conductivity = np.gradient(list(map(lambda x: Current(x),temps2)), temps2)
    ax[4].plot(list(map(lambda x: Volt(x) / 10 ** (-3),temps2)), conductivity, label = "set "+str(num))

    ax[4].set_title("Проводимость от напряжения")
    ax[4].set_xlabel("Volt(Temp)")
    ax[4].set_ylabel("Conductivity")
    ax[4].ticklabel_format(style = 'scientific', scilimits = [-1, 1])
    ax[4].legend(loc = 4)
    # Мощность необходимая для перехода линий в нормальное состояние

    def TempInt(temp):
        return temp * Int(temp)

    def CurrCool(T_bath, T_Crit):
        return (gamma * delta(T_bath) / e_0 + np.sqrt( (gamma * delta(T_bath) / e_0) ** 2 + 64 * k_B ** 2 / e_0 ** 2 * R_g / R_lead * quad(TempInt, T_bath, T_Crit)[0])) / (2 * R_g)

    t_baths = np.arange(0.0001, 2, 0.001)
    tbtc = list(map(lambda t: t / T_C, t_baths))
    ps = list(map(lambda tb: CurrCool(tb, T_C) / 10 ** (-6), t_baths))
    ax[5].plot(tbtc, ps, label = "set "+str(num))

    ax[5].set_title("Мощность необходимая для перехода линий в нормальное состояние")
    ax[5].set_xlabel(r'$\frac{T_b}{T_c}$')
    ax[5].set_ylabel("P")
    ax[5].ticklabel_format(style = 'scientific')
    ax[5].legend(loc = 4)

nums = [i for i in range(nsets)]

list(map(lambda x, y: DrawGraphs(*x, y), vars, nums))

plt.show()
