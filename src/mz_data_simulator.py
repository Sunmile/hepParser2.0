import time

import numpy as np
import pandas as pd
import pickle as pk
import os
from src.enumerate import get_comp_pk
import matplotlib.pyplot as plt
from brainpy import isotopic_variants
from src.utils import save_pk
# from scipy.stats import poisson
from scipy.stats import gamma
from scipy.stats import norm
import random
import settings


def get_all_hp_component():
    min_dp = settings.min_dp
    max_dp = settings.max_dp
    the_HP = get_comp_pk('./../data/enox_', min_dp, max_dp)
    mass_list = the_HP[4]
    dict_mass_atom = the_HP[1]
    return mass_list, dict_mass_atom


def get_iso_distribution(mass_list, dict_mass_atom):
    iso_num = settings.iso_num
    result = []
    all_mass = sorted(mass_list)
    flag = all_mass[0]
    for x in all_mass:
        if x - flag > 10:
            tmp_com = dict_mass_atom[x][0]
            tmp_lwh = {}
            tmp_lwh['H'] = tmp_com[0]
            tmp_lwh['C'] = tmp_com[1]
            tmp_lwh['N'] = tmp_com[2]
            tmp_lwh['O'] = tmp_com[3]
            tmp_lwh['S'] = tmp_com[4]
            tmp_distribution = isotopic_variants(tmp_lwh, npeaks=iso_num, charge=0)
            tmp_point = [x]
            for peak in tmp_distribution:
                tmp_point.append(peak.intensity)
            result.append(tmp_point)
            flag = x
    return result


def fit_all_point(point_list, degree=4):
    iso_num = settings.iso_num
    x = []
    y_list = []
    para_list = []
    fit_y_list = []
    for p in point_list:
        x.append(p[0])
    for i in range(1, iso_num + 1):
        tmp_y = []
        for p in point_list:
            tmp_y.append(p[i])
        y_list.append(tmp_y)
        tmp_z = np.polyfit(x, tmp_y, degree)
        tmp_parm = np.poly1d(tmp_z)
        tmp_y_fit = tmp_parm(x)
        para_list.append(tmp_parm)
        fit_y_list.append(tmp_y_fit)
    plt.figure(figsize=(12, 10), dpi=200)
    color_list = ['#f4a587', '#70b543', '#4bbcec', '#fdb514', '#ea5f9a', '#c26567', '#9f9d58',
                  '#9ba439', '#57b67a', '#979998', '#51632d', '#4b5066', '#98b1dc', '#e1a626']
    for i in range(len(y_list)):
        tmp_y = y_list[i]
        tmp_y_fit = fit_y_list[i]
        plt.scatter(x, tmp_y, color=color_list[i], s=5)
        plt.plot(x, tmp_y_fit, color=color_list[i], label=str(i + 1))
    plt.legend()
    plt.savefig('./../data/fit.png')
    plt.show()
    return para_list


def get_fit_iso_pk(dir, degree=4):
    min_dp = settings.min_dp
    max_dp = settings.max_dp
    iso_num = settings.iso_num
    dir = dir + '_' + str(min_dp) + '_' + str(max_dp) + '_' + str(iso_num) + '_' + str(degree) + '.pk'
    if os.path.exists(dir):
        with open(dir, 'rb') as f:
            data = pk.load(f)
        return data
    else:
        mass_list, dict_mass_atom = get_all_hp_component()
        point_list = get_iso_distribution(mass_list, dict_mass_atom)
        para_list = fit_all_point(point_list, degree)
        save_pk(para_list, dir)
        return para_list


def get_noise_rand_list(col_num):
    res_list = []
    for i in range(col_num):
        x = random.uniform(0, 1)
        res_list.append(x * settings.noise_int_ratio)
    return res_list


def get_range_gamma(a, b):
    for i in range(int(a * b), 100):
        tmp_p = gamma.pdf(i, a=a, scale=b)
        if tmp_p < 0.05:
            return i
    return 100


def get_rand_list(mu, std, num, x_left_bias, x_right_bias, y_bias, ratio_list, noise_mz_num, dis='norm'):
    noise_mz_num_std = settings.noise_mz_num_std
    result_list = []
    # flag_list = []  # 与result_list对应， 用于记录 其中元素是有效峰(1)还是噪音峰(0)
    count_list = []  # 记录每两个同位素峰之间隔了多少个噪音峰
    base_list = []
    if dis == 'norm' or dis == 'noise_norm' or dis == 'norm_overlap':
        tmp_inf = mu - (3 + x_left_bias) * std
        tmp_sup = mu + (3 + x_right_bias) * std
        interval = np.round((tmp_sup - tmp_inf) / num, 5)
        print('interval:', interval)
        max_p = norm.pdf(mu, loc=mu, scale=std)
        x = tmp_inf
        for i in range(num):
            p = norm.pdf(x, loc=mu, scale=std)
            base_list.append(p / max_p)
            x = x + interval
    elif dis == 'gamma':
        alpha = random.uniform(1.5, 10)
        beta = random.uniform(0.2, 5)
        tmp_inf = 0 - x_left_bias
        tmp_sup = get_range_gamma(alpha, beta) + x_right_bias
        if tmp_sup - tmp_inf < 1:
            tmp_inf = tmp_inf - (np.abs(tmp_sup - tmp_inf) + 1) / 2
            tmp_sup = tmp_sup + (np.abs(tmp_sup - tmp_inf) + 1) / 2
        interval = np.round((tmp_sup - tmp_inf) / num, 5)
        print('interval:', interval)
        x = tmp_inf
        max_p = 0
        tmp_list = []
        for i in range(num):
            p = gamma.pdf(x, a=alpha, scale=std)
            if p > max_p:
                max_p = p
            tmp_list.append(p)
            x = x + interval
        for i in range(num):
            base_list.append(tmp_list[i] / max_p)
    elif dis == 'noise':
        base_list = np.random.uniform(0, 1, num)
    elif dis == 'noise_rand':
        base_list = np.random.uniform(0, 1, num)
        ratio_list = np.random.uniform(0,1,len(ratio_list))
    else:
        print('Now only support "norm","norm_overlap", "gamma", "noise", "noise_rand" and "noise_norm" distribution')
        return
    if dis == 'noise_norm':
        random.shuffle(ratio_list)
    if dis == 'norm_overlap':
        start_1 = random.randint(0, len(ratio_list)-1)
        start_2 = random.randint(0, len(ratio_list)-1)
        step = random.randint(1, 5)
        choose = random.choice([0, 1])
        ratio_list_1 = ratio_list
        ratio_list_2 = ratio_list
        if choose == 0:
            for i in range(start_1, len(ratio_list), step):
                ratio_list_1[i] += ratio_list_2[start_2]
                start_2 += 1
                if start_2 >=len(ratio_list):
                    break
        else:
            for i in range(start_2, len(ratio_list), step):
                ratio_list_1[start_1] += ratio_list_2[i]
                start_1 += 1
                if start_1 >=len(ratio_list):
                    break
        max_ratio = max(ratio_list_1)
        if max_ratio <= 0:
            max_ratio = 1
        ratio_list = [x / max_ratio for x in ratio_list_1]
    if len(base_list) != num:
        print('no match columns', num, len(base_list))
    save_base_list = base_list
    for charge in range(1, settings.max_charge + 1):
        tmp_noise_num = max(1, int(noise_mz_num * (1 / charge)))
        tmp_res_list = []
        # tmp_flag_list = []
        tmp_count_list = []
        count_norm = 0
        for i in range(len(ratio_list)):
            if dis == 'noise_norm':
                x = random.randint(0, 1)
                if x == 1 and count_norm < settings.max_num_norm:
                    base_list = save_base_list
                    count_norm += 1
                else:
                    base_list = np.random.uniform(0, 1, num)
            one_list = []
            for j in range(len(base_list)):
                tmp_bias = random.uniform(-y_bias,y_bias)
                new_int = (base_list[j] + tmp_bias) * ratio_list[i]
                if new_int < 0:
                    new_int = 0
                one_list.append(new_int)
            one_list = random_set_0(one_list, settings.p_set0)
            tmp_res_list.append(one_list)
            # tmp_flag_list.append(1)
            if i < len(ratio_list) - 1:
                noise_inf = int(tmp_noise_num * (1 - noise_mz_num_std))
                noise_sup = int(tmp_noise_num * (1 + noise_mz_num_std))
                noise_num = random.randint(noise_inf, noise_sup)
                for k in range(noise_num):
                    tmp_list = get_noise_rand_list(num)
                    tmp_list = random_set_0(tmp_list, settings.p_set0)
                    tmp_res_list.append(tmp_list)
                    # tmp_flag_list.append(0)
                tmp_count_list.append(noise_num)
        result_list.append(tmp_res_list)
        # flag_list.append(tmp_flag_list)
        count_list.append(tmp_count_list)
    return result_list, count_list


def get_number_by_pro(number_list, pro_list):
    """
    :param number_list:数字列表
    :param pro_list:数字对应的概率列表
    :return:按概率从数字列表中抽取的数字
    """
    # 用均匀分布中的样本值来模拟概率
    x = random.uniform(0, 1)
    # 累积概率
    cum_pro = 0.0
    # 将可迭代对象打包成元组列表
    for number, number_pro in zip(number_list, pro_list):
        cum_pro += number_pro
        if x < cum_pro:
            return number


# 按正态分布选择元素
def normal_choice(lst, mean=None, stddev=None):
    if mean is None:
        # if mean is not specified, use center of list
        mean = (len(lst) - 1) / 2

    if stddev is None:
        # if stddev is not specified, let list be -3 .. +3 standard deviations
        stddev = len(lst) / 6

    while True:
        index = int(random.normalvariate(mean, stddev) + 0.5)
        if 0 <= index < len(lst):
            return lst[index]


def random_set_0(value_list, p=0.1):
    tmp_list = [0, 1]
    tmp_p = [p, 1 - p]
    new_list = []
    for i in range(len(value_list)):
        num = get_number_by_pro(tmp_list, tmp_p)
        new_list.append(num * value_list[i])
    return new_list


def random_mismatch(value_list, num=11):
    column = len(value_list)
    new_arr = np.zeros([num, column])
    row_list = list(range(num))
    for i in range(column):
        ind = normal_choice(row_list)
        new_arr[ind, i] = value_list[i]
    return new_arr


def get_iso_ratio(mz, para_list):
    ratio_list = []
    for i in range(min(len(para_list), settings.iso_num)):
        tmp_parm = para_list[i]
        ratio_list.append(tmp_parm(mz))
    max_ratio = max(ratio_list)
    ratio_list = [(one / max_ratio) for one in ratio_list]
    return ratio_list


def get_main_mzlist(tmp_mass, count_list, charge):
    one_mz_list = []
    mz = (tmp_mass - charge) / charge
    one_mz_list.append(mz)
    for j in range(len(count_list)):
        tmp_mz_list = []
        mz_0 = mz + j / charge
        for k in range(count_list[j]):
            x = random.uniform(0, 1 / charge)
            tmp_mz = mz_0 + x
            tmp_mz_list.append(tmp_mz)
        one_mz_list.extend(sorted(tmp_mz_list))
        one_mz_list.append(mz_0 + 1 / charge)
    return one_mz_list


def get_random_mzlist(tmp_mass, row_num_list, count_list, charge):
    all_mz_list = []
    main_mz_list = get_main_mzlist(tmp_mass, count_list, charge)
    mz = (tmp_mass - charge) / charge
    if len(main_mz_list) != len(row_num_list):
        print('unmatched mz num and distribution num:', len(main_mz_list), len(row_num_list))
        print('count_list:', count_list)
        print('sum of count:', sum(count_list))
    win = mz * settings.ppm * 0.000001
    for j in range(len(row_num_list)):
        one_row_num = row_num_list[j]
        mz_0 = main_mz_list[j]
        for k in range(one_row_num):
            tmp_v = norm.rvs(loc=0, scale=0.2, size=1)[0]
            if tmp_v > 1:
                tmp_v = 1
            if tmp_v < -1:
                tmp_v = -1
            mz_shift = tmp_v * win
            new_mz = mz_0 + mz_shift
            all_mz_list.append(np.round(new_mz, 5))
    return all_mz_list


def generate_positive_data(dir, para_list, dis='norm'):
    min_mass = settings.min_mass
    max_mass = settings.max_mass
    dir = dir + '/{dis}_{i}_{j}.csv'
    all_mass_list = list(range(min_mass, max_mass))
    mu_list = list(range(10))
    std_list = list(np.arange(1, 5, 0.1))
    x_bias_list = list(np.arange(-1, 3, 0.1))
    y_bias_list = list(np.arange(0, 0.5, 0.001))
    sample_num = settings.hp_num
    max_row_num = settings.max_row_num
    max_columns_num = settings.max_columns_num
    for i in range(sample_num):
        sample_result = []
        tmp_mz = random.choice(all_mass_list)
        ratio_list = get_iso_ratio(tmp_mz, para_list)
        tmp_mu = random.choice(mu_list)
        tmp_std = random.choice(std_list)
        tmp_xleft_bias = random.choice(x_bias_list)
        tmp_xright_bias = random.choice(x_bias_list)
        tmp_y_bias = random.choice(y_bias_list)
        columns_num = random.choice(range(20, max_columns_num))
        noise_mz_num = random.choice(range(10, settings.max_noise_mz_num))
        result_list, count_list = \
            get_rand_list(mu=tmp_mu,
                          std=tmp_std,
                          num=columns_num,
                          x_left_bias=tmp_xleft_bias,
                          x_right_bias=tmp_xright_bias,
                          y_bias=tmp_y_bias,
                          ratio_list=ratio_list,
                          noise_mz_num=noise_mz_num,
                          dis=dis
                          )
        for k in range(settings.max_charge):
            tmp_result_list = result_list[k]
            tmp_count_list = count_list[k]
            row_num_list = []
            for j in range(len(tmp_result_list)):
                row_num = random.choice(range(1, max_row_num))
                row_num_list.append(row_num)
                new_result = random_mismatch(tmp_result_list[j], row_num)
                if j == 0:
                    sample_result = new_result
                else:
                    sample_result = np.concatenate([sample_result, new_result], axis=0)
            if dis == 'noise_rand':
                sel = random.choice([0,1])
                if sel == 0:
                    random.shuffle(sample_result)
            all_mzlist = get_random_mzlist(tmp_mz, row_num_list, tmp_count_list, charge=k + 1)
            mz_array = np.array(all_mzlist).reshape([-1, 1])
            one_sample = np.concatenate([mz_array, sample_result], axis=1)
            one_sp_df = pd.DataFrame(one_sample)
            one_sp_df = one_sp_df.sort_values(by=0)
            new_dir = dir.format(dis=dis, i=i, j=k + 1)
            one_sp_df.to_csv(new_dir, index=None, header=None)
        print('Has generate data num:', i)


if __name__ == '__main__':
    dir_fit = './../data/fit'
    dir_pos = './../simulated_data/positive'
    para_list = get_fit_iso_pk(dir_fit)
    start = time.time()
    # Now only support "norm","norm_overlap", "gamma", "noise","noise_rand" and "noise_norm" distribution
    generate_positive_data(dir_pos, para_list, dis='norm')
    print('Time cost', time.time() - start)

    # color_list = ['#f4a587', '#70b543', '#4bbcec', '#fdb514', '#ea5f9a', '#c26567', '#9f9d58',
    #               '#9ba439', '#57b67a', '#979998', '#51632d', '#4b5066', '#98b1dc', '#e1a626']
    # #
    # p_list = []
    # # p3_list = []
    # #
    # alpha = 10
    # beta = 5
    # #
    # x_list = np.arange(0,100,0.1)
    # for x in x_list:
    #     p = gamma.pdf(x, a=alpha,scale=beta)
    #     # p3 = norm.pdf(x, loc=20,scale=4)
    #     p_list.append(p)
    #     # p3_list.append(p3)
    # plt.figure(figsize=(6,4),dpi=200)
    # plt.plot(x_list,p_list,c=color_list[2],label='1,1')
    # plt.legend()
    # plt.show()
