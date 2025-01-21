import numpy as np


def combo(x: list, y1=None, y2=None):
    """
    Функция объениняет n и k в один масив, чтобы в дальнейшем использовать их для одновременного фитирования
    :param x: mg value like of list
    :param y1: n value like of list
    :param y2: k value like of list
    :return: 2 list if the values are set mg, n and k, or 1 list if only the value of mg is set
    """
    combo_x = np.append(x, x)
    combo_y = np.append(y1, y2)
    if (y1 is None) and (y2 is None):
        return combo_x
    else:
        return combo_x, combo_y


def pw_func_n(mg, mg1, nm, nb, nu, km, kb, ku):
    return np.piecewise(mg, [mg <= mg1, mg > mg1],
                        [
                            lambda mg: nm + nb * mg,
                            lambda mg: nm + nb * mg1 + nu * (mg - mg1),
                        ])


def pw_func_k(mg, mg1, nm, nb, nu, km, kb, ku):
    return np.piecewise(mg, [mg <= mg1, mg > mg1],
                        [
                            lambda mg: km + kb * mg,
                            lambda mg: km + kb * mg1 + ku * (mg - mg1),
                        ])


def pw_combo(mg, mg1, nm, nb, nu, km, kb, ku):
    combo_data = combo(mg)
    data1 = combo_data[:len(mg)]
    data2 = combo_data[len(mg):]

    res1 = pw_func_n(data1, mg1, nm, nb, nu, km, kb, ku)
    res2 = pw_func_k(data2, mg1, nm, nb, nu, km, kb, ku)

    return np.append(res1, res2)
