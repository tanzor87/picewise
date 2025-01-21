import numpy as np
from pwlf import PiecewiseLinFit
from scipy import optimize
import pwlf
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path


def nearest(lst, target):
    return min(lst, key=lambda x: abs(x - target))


def get_nk_single_t_from_excel(file_name: Path, temperature):
    nk_file = pd.ExcelFile(file_name)
    data_nk_single_t = nk_file.parse(sheet_name=temperature, header=[0, 1], index_col=0)
    return data_nk_single_t


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


def initial_parameters(data, temperature, num_breaks, target_frreq):
    frequencies = data.index.to_numpy()
    frequency = nearest(frequencies, freq_target)

    df_n_single_f = data.loc[frequency, ['n_red']].dropna()

    n_single_f = df_n_single_f.loc['n_red'].to_numpy()
    k_single_f = data.loc[frequency, ['k_red']].dropna().to_numpy()
    moisture = df_n_single_f.loc['n_red'].index.to_numpy()

    pwlf1 = pwlf.PiecewiseLinFit(moisture, n_single_f)
    pwlf1.use_custom_opt(num_breaks + 1)
    pwlf2 = pwlf.PiecewiseLinFit(moisture, k_single_f)
    pwlf2.use_custom_opt(num_breaks + 1)


if __name__ == '__main__':
    file_name = Path(r'D:\Python\picewise\picewise\Data\nk_red_rev_p18_nk_mol.xlsx')
    temperature = '-30'
    num_breaks = 2
    freq_target = 0.7e9

    data = get_nk_single_t_from_excel(file_name, temperature)

    # moisture = data['n_red'].columns.to_numpy()
    frequencies = data.index.to_numpy()
    frequency = nearest(frequencies, freq_target)

    df_n_single_f = data.loc[frequency, ['n_red']].dropna()
    n_single_f = df_n_single_f.loc['n_red'].to_numpy()
    k_single_f = data.loc[frequency, ['k_red']].dropna().to_numpy()

    moisture = df_n_single_f.loc['n_red'].index.to_numpy()

    pwlf1 = pwlf.PiecewiseLinFit(moisture, n_single_f)
    pwlf1.fit(num_breaks + 1)
    pwlf2 = pwlf.PiecewiseLinFit(moisture, k_single_f)
    pwlf2.fit(num_breaks + 1)

    mg_by_n = pwlf1.fit_breaks[1: -1]
    mg_by_k = pwlf2.fit_breaks[1: -1]
    print(f'mg_by_n = {mg_by_n}')
    print(f'mg_by_k = {mg_by_k}')
    mg_mean = (mg_by_n + mg_by_k) / 2
    print(f'mg_mean = {mg_mean}')

    moisture_combo, nk_single_f = combo(moisture, n_single_f, k_single_f)

    my_combo_pwlf = pwlf.PiecewiseLinFit(moisture_combo, nk_single_f)

    my_combo_pwlf.fit(num_breaks + 1)

    print(f'mg_combo = {my_combo_pwlf.fit_breaks[1: -1]}')

    # my_pwlf_n = PiecewiseLinFit(moisture, )

    #     ===================================================

    pwlf1 = pwlf.PiecewiseLinFit(moisture, n_single_f)
    pwlf1.use_custom_opt(num_breaks + 1)
    pwlf2 = pwlf.PiecewiseLinFit(moisture, k_single_f)
    pwlf2.use_custom_opt(num_breaks + 1)


    def mse(breakpoint):
        # print(breakpoint, breakpoint.shape, breakpoint.ndim)
        ssr1 = pwlf1.fit_with_breaks_opt(breakpoint)
        mse1 = ssr1 / pwlf1.n_data
        ssr2 = pwlf2.fit_with_breaks_opt(breakpoint)
        mse2 = ssr2 / pwlf2.n_data
        return mse1 + mse2


    mg1_2 = np.zeros(2)
    # mg1_2[0] = max(moisture) / 3
    mg1_2[0] = min(moisture) + (max(moisture) - min(moisture)) / 3
    # mg1_2[1] = max(moisture) * 2 / 3
    mg1_2[1] = min(moisture) + 2 * (max(moisture) - min(moisture)) / 3
    print(mg1_2)
    bounds = np.array([[moisture.min(), moisture.max()]])
    breakpoints = optimize.minimize(mse, mg1_2, bounds=bounds, method='L-BFGS-B')

    print(breakpoints['x'])
