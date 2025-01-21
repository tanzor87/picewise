from pathlib import Path
import functions as fnc
import matplotlib.pyplot as plt
import pandas as pd
from scipy import optimize
from pw_3_water_nk import pw_func_n, pw_func_k, pw_combo
import numpy as np


np.set_printoptions(precision=4)

path_file_name = Path(r'D:\Python\picewise\picewise\Data\nk_red_rev_p18_nk_mol.xlsx')
file_name = path_file_name.stem
temperature = '-20'
num_breaks = 2
freq_target = 1.1e9
fixed = 0              # 0 - без фиксации mg1 и mg2, 1 - с фиксацией
n = 0                  # Количество влажностей с конца, которые будут удалены, если это необходимо

# Считываем Excel файл при нужной температуре
data = fnc.get_nk_single_t_from_excel(path_file_name, temperature)
# Получаем зависимость от влажности на интересующей нас частоте
target_data = fnc.chose_data_target_freq(data, freq_target, last_del=n)
# Получаем первоначальные значения для точек излома
initial_breaks = fnc.initial_mg_ts(target_data, num_breaks)
# Получаем остальные необходимые параметры модели
n_intercept, k_intercept, n_slope, k_slope, predict_moisture, n_fit, k_fit = fnc.initial_slopes(target_data,
                                                                                                initial_breaks)
# initial = initial_breaks + n_intercept[0] + n_slope + k_intercept[0] + k_slope
initial = np.concatenate([initial_breaks, [n_intercept[0]], n_slope, [0.0], k_slope])

print('Первоначальные параметры для фитирования: ')
print(
    f'mg_ts = {initial_breaks},\n'
    f'n_slopes = {n_slope},\n'
    f'k_slopes = {k_slope},\n'
    f'n_intercept = {n_intercept[0]},\n'
    f'k_intercept = {k_intercept[0]}'
)

# Список частот
frequencies = fnc.get_frequency(data)
freq_target = fnc.nearest(frequencies, freq_target)
# Список температур
temperatures = fnc.get_temperatures(path_file_name)

# Создаем пустой датафрейм для полученных параметров фита
columns = pd.MultiIndex.from_product([['mg1', 'mg2', 'nm', 'nb', 'nt', 'nu', 'km', 'kb', 'kt', 'ku'],
                                      ['value', 'error']], names=['parametr', 'v/e'])
rows = pd.Index(frequencies, name='Frequency')
fit_params = pd.DataFrame(index=rows, columns=columns)
# print(fit_params)

# name_create = "1.xlsx"
if fixed == 0:
    # ----цикл без фиксации mg1 и mg1-----------------------------------------------------------------
    for frequency in frequencies:
        moisture, n_single_f, k_single_f = fnc.chose_data_target_freq(data, frequency, last_del=n)

        fit_param, e = optimize.curve_fit(
            pw_combo,
            moisture,
            fnc.combo(
                moisture,
                n_single_f,
                k_single_f)[1],
                initial,
                bounds=((0, 0, -np.inf, -np.inf, -np.inf,
                         -np.inf, 0, 0.03, -np.inf, -np.inf),
                        (1, 1, np.inf, np.inf, np.inf,
                         np.inf, np.inf, np.inf, np.inf, np.inf)))
        p_err = np.sqrt(np.diag(e))
        fit_params.loc[frequency, (slice(None), 'value')] = fit_param
        fit_params.loc[frequency, (slice(None), 'error')] = p_err
        # print(fit_param)
        # initial = fit_param

elif fixed == 1:
    # ----цикл с фиксацией mg1 и mg1-----------------------------------------------------------------------------------
    eps = 0.00001
    for frequency in frequencies:
        moisture, n_single_f, k_single_f = fnc.chose_data_target_freq(data, frequency, last_del=n)

        fit_param, e = optimize.curve_fit(
            pw_combo,
            moisture,
            fnc.combo(
                moisture,
                n_single_f,
                k_single_f)[1],
                initial,
                bounds=((initial[0] - eps, initial[1] - eps, -np.inf, -np.inf, -np.inf,
                         -np.inf, 0, 0.03, -np.inf, -np.inf),
                        (initial[0] + eps, initial[1] + eps, np.inf, np.inf, np.inf,
                         np.inf, np.inf, np.inf, np.inf, np.inf)))
        p_err = np.sqrt(np.diag(e))
        fit_params.loc[frequency, (slice(None), 'value')] = fit_param
        fit_params.loc[frequency, (slice(None), 'error')] = p_err
        # initial = fit_param

name_create = Path(f'Result/spectra_of_params_{file_name}.xlsx')

# #  ------- записываем спектры в exel файл --------
# if Path.is_file(name_create):
#     with pd.ExcelWriter(name_create, engine='openpyxl', mode='a') as writer:
#         fit_params.astype(float).round(6).to_excel(writer, sheet_name=temperature)
# else:
#     with pd.ExcelWriter(name_create, engine='openpyxl') as writer:
#         fit_params.astype(float).round(6).to_excel(writer, sheet_name=temperature)


'''Извлекаем данные для построения графиков'''

mg1 = fit_params.loc[frequencies, (['mg1'], 'value')].T.to_numpy()[0]
mg2 = fit_params.loc[frequencies, (['mg2'], 'value')].T.to_numpy()[0]

nb = fit_params.loc[frequencies, (['nb'], 'value')].T.to_numpy()[0]
nt = fit_params.loc[frequencies, (['nt'], 'value')].T.to_numpy()[0]
nu = fit_params.loc[frequencies, (['nu'], 'value')].T.to_numpy()[0]

kb = fit_params.loc[frequencies, (['kb'], 'value')].T.to_numpy()[0]
kt = fit_params.loc[frequencies, (['kt'], 'value')].T.to_numpy()[0]
ku = fit_params.loc[frequencies, (['ku'], 'value')].T.to_numpy()[0]


mg1_err = fit_params.loc[frequencies, (['mg1'], 'error')].T.to_numpy()[0]
mg2_err = fit_params.loc[frequencies, (['mg2'], 'error')].T.to_numpy()[0]
nb_err = fit_params.loc[frequencies, (['nb'], 'error')].T.to_numpy()[0]
nt_err = fit_params.loc[frequencies, (['nt'], 'error')].T.to_numpy()[0]
nu_err = fit_params.loc[frequencies, (['nu'], 'error')].T.to_numpy()[0]
kb_err = fit_params.loc[frequencies, (['kb'], 'error')].T.to_numpy()[0]
kt_err = fit_params.loc[frequencies, (['kt'], 'error')].T.to_numpy()[0]
ku_err = fit_params.loc[frequencies, (['ku'], 'error')].T.to_numpy()[0]

param_single_f = fit_params.loc[freq_target, (slice(None), 'value')].to_numpy()
param_single_f = [round(float(i), 4) for i in param_single_f]
print(f'Найденные параметры на частоте f = {round(freq_target / 1e9, 3)} ГГц:\n{param_single_f}')

# # -------- модельные значения n и k --------
moisture = target_data[0]
n_single_f = target_data[1]
k_single_f = target_data[2]

mg_mod = np.linspace(0, max(moisture), 200)
n_mod = pw_func_n(mg_mod, *param_single_f)
k_mod = pw_func_k(mg_mod, *param_single_f)

# Графики
fig1, ([[ax11, ax12], [ax21, ax22]]) = plt.subplots(2, 2)
fig1.set_size_inches((15, 10))
plt.gcf().subplots_adjust(left=0.06, bottom=0.06, right=0.95, top=0.97, wspace=0.177, hspace=0.18)

ax11.errorbar(frequencies, mg1, yerr=mg1_err, fmt='o')
ax11.errorbar(frequencies, mg2, yerr=mg2_err, fmt='s')

ax11.set_title(f'T = {temperature} °C')
ax11.set_xlabel('Частота, Гц')
ax11.set_ylabel('mg1, mg2')
ax11.set_xscale('log')
ax11.legend(['mg1', 'mg2'])
ax11.grid()

ax12.plot(frequencies, nb)
ax12.plot(frequencies, nt)
ax12.plot(frequencies, nu)

# ax12.errorbar(frequencies, nb, yerr=nb_err, fmt='o')
# ax12.errorbar(frequencies, nt, yerr=nt_err, fmt='o')
# ax12.errorbar(frequencies, nu, yerr=nu_err, fmt='o')

ax12.set_title(f'T = {temperature} °C')
ax12.set_xlabel('Частота, Гц')
ax12.set_ylabel('Нормированное n категорий воды')
# plt.text(0.3e9, 18.0, f'T = {temperature} °C')
ax12.set_xscale('log')
ax12.legend(['nb', 'nt', 'nu'])
ax12.grid()

ax22.plot(frequencies, kb)
ax22.plot(frequencies, kt)
ax22.plot(frequencies, ku)

ax22.set_title(f'T = {temperature} °C')
ax22.set_xlabel('Частота, Гц')
ax22.set_ylabel(r'Нормированное $\kappa$ категорий воды')
# plt.text(0.3e9, 14.0, f'T = {temperature} °C')
ax22.set_xscale('log')
ax22.legend(['kb', 'kt', 'ku'])
ax22.grid()

ax21.plot(moisture, n_single_f, 'ob', label='(n-1) / ro')
ax21.plot(mg_mod, n_mod, '-b')
ax21.plot(moisture, k_single_f, 'og', label='k / ro')
ax21.plot(mg_mod, k_mod, '-g')
ax21.set_xlabel('Влажность, mg')
ax21.set_ylabel('(n-1) / ro,  k / ro')
ax21.set_title(f'Частота - {round(freq_target / 1e9, 3)} ГГц, T = {temperature} °C')
# plt.text(0.1, 0.9, f'T = {temperature} °C')
ax21.legend()
ax21.grid()

plt.show()
