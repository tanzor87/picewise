import numpy as np
from scipy import optimize
import pwlf
import pandas as pd
from pathlib import Path


def nearest(lst, target):
    return min(lst, key=lambda x: abs(x - target))


def get_nk_single_t_from_excel(file_name: Path, temperature):
    nk_file = pd.ExcelFile(file_name)
    data_nk_single_t = nk_file.parse(sheet_name=temperature, header=[0, 1], index_col=0)
    return data_nk_single_t


def get_temperatures(file_name: Path):
    data_file = pd.ExcelFile(file_name)
    temperatures = [sheet for sheet in data_file.sheet_names]
    return temperatures


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


def chose_data_target_freq(data, target_freq, last_del=0):
    frequencies = data.index.to_numpy()
    frequency = nearest(frequencies, target_freq)

    df_n_single_f = data.loc[frequency, ['n_red']].dropna()

    n_single_f = df_n_single_f.loc['n_red'].to_numpy()
    k_single_f = data.loc[frequency, ['k_red']].dropna().to_numpy()
    moisture = df_n_single_f.loc['n_red'].index.to_numpy()

    if last_del != 0:
        n_single_f = n_single_f[:-last_del]
        k_single_f = k_single_f[:-last_del]
        moisture = moisture[:-last_del]

    return moisture, n_single_f, k_single_f


def get_frequency(data):
    frequencies = data.index.to_numpy()
    return frequencies


def initial_mg_ts(data, num_breaks):
    moisture = data[0]
    n_single_f = data[1]
    k_single_f = data[2]

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

    mg_ts = np.zeros(2)
    mg_ts[0] = min(moisture) + (max(moisture) - min(moisture)) / 3
    mg_ts[1] = min(moisture) + 2 * (max(moisture) - min(moisture)) / 3

    bounds = np.array([[moisture.min(), moisture.max()]])
    breakpoints = optimize.minimize(mse, mg_ts, bounds=bounds, method='L-BFGS-B')
    breakpoints = breakpoints['x']

    return breakpoints


def initial_slopes(data, mg_ts, num_points = 400):
    moisture = data[0]
    n_single_f = data[1]
    k_single_f = data[2]

    pwlf_n = pwlf.PiecewiseLinFit(moisture, n_single_f)
    pwlf_k = pwlf.PiecewiseLinFit(moisture, k_single_f)

    breakpoints = [min(moisture), mg_ts[0], mg_ts[1], max(moisture)]

    pwlf_n.fit_with_breaks(breakpoints)
    pwlf_k.fit_with_breaks(breakpoints)

    n_slopes = pwlf_n.slopes
    k_slopes = pwlf_k.slopes

    # breaks = pwlf_n.fit_breaks

    n_intercepts = pwlf_n.intercepts
    k_intercepts = pwlf_k.intercepts

    x_hat = np.linspace(min(moisture), max(moisture), num_points)
    n_fit = pwlf_n.predict(x_hat)
    k_fit = pwlf_k.predict(x_hat)

    return n_intercepts, k_intercepts, n_slopes, k_slopes, x_hat, n_fit, k_fit


# if __name__ == '__main__':
    # file_name = Path(r'D:\Python\picewise\picewise\Data\nk_red_rev_p18_nk_mol.xlsx')
    # temperature = '-30'
    # num_breaks = 2
    # freq_target = 0.7e9
    #
    # data = get_nk_single_t_from_excel(file_name, temperature)
    #
    # target_data = chose_data_target_freq(data, freq_target)
    #
    # initial_breaks = initial_mg_ts(target_data, num_breaks)
    #
    # print(initial_breaks)
    #
    # n_intercept, k_intercept, n_spope, k_slope = initial_slopes(target_data, initial_breaks)
    #
    # print(
    #     f'n_slopes = {n_spope},\n'
    #     f'k_slopes = {k_slope},\n'
    #     f'n_intercept = {n_intercept[0]},\n'
    #     f'k_intercept = {k_intercept[0]}'
    #     )

