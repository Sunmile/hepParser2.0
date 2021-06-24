import numpy as np
import pandas as pd
import copy
from src.enumerate import *
import pickle as pk


def read_spectra(dir):
    MZ_int = []
    max_int = 0
    with open(dir, 'r') as f:
        line = f.readline()
        for i in range(5):
            line = f.readline()
        while not 'END' in line:
            info = line.strip('\n').split(' ')
            MZ = float(info[0])
            intensity = float(info[1])
            if MZ >= 0:  # 350
                if intensity > max_int:
                    max_int = intensity
                MZ_int.append([np.round(MZ, 5), intensity])
            line = f.readline()

        for i in range(len(MZ_int)):
            reletive_int = MZ_int[i][1] / max_int
            MZ_int[i].append(np.round(reletive_int, 6))
        print('Max Intensity:' + str(max_int))
    return MZ_int, max_int


def get_sp_input(Mz_list):
    MZ_int = []
    max_int = 0
    for i in Mz_list:
        MZ = float(i[0])
        intensity = float(i[1])
        if MZ >= 0:  # 350
            if intensity > max_int:
                max_int = intensity
            MZ_int.append([MZ, intensity])
    for i in range(len(MZ_int)):
        reletive_int = MZ_int[i][1] / max_int
        MZ_int[i].append(np.round(reletive_int, 6))
    return MZ_int, max_int


def filter_component(all_list, lowth, highth):
    result_list = []
    for x in all_list:
        mass, tmp = get_component_mass(x)
        if lowth <= mass <= highth:
            result_list.append(x)
    return result_list


def save_file(data, dir):
    f = open(dir, 'w', newline='')
    for i in range(len(data)):
        f.writelines(str(data[i][0]) + ' ' + str(data[i][1]))
        f.writelines('\n')
    f.close()


def save_pk(data, dir):
    with open(dir, 'wb') as f:
        pk.dump(data, f)


def cos_sim(vector_a, vector_b):
    """
    计算两个向量之间的余弦相似度
    :param vector_a: 向量 a
    :param vector_b: 向量 b
    :return: sim [0,1]
    :return: cos [-1,1] 由于本任务的特殊性，cos [0, 1]
    """
    vector_a = np.mat(vector_a)
    vector_b = np.mat(vector_b)
    num = float(vector_a * vector_b.T)
    denom = np.linalg.norm(vector_a) * np.linalg.norm(vector_b)
    cos = num / denom
    # sim = 0.5 + 0.5 * cos
    return cos


def get_KL(v1, v2, is_P=0):
    if is_P == 0:
        v1 = np.array(v1)
        v2 = np.array(v2)
        p_v1 = v1 / sum(v1)
        p_v2 = v2 / sum(v2)
    else:
        p_v1 = v1
        p_v2 = v2
    KL = 0
    for i in range(len(v1)):
        if p_v1[i] > 0 and p_v2[i] > 0:
            KL += p_v1[i] * np.log2(p_v1[i] / p_v2[i])
    return KL


def get_JS(v1, v2, w1=0.5, w2=0.5):
    v1 = np.array(v1)
    v2 = np.array(v2)

    p_v1 = v1 / sum(v1)
    p_v2 = v2 / sum(v2)
    p_avg = (p_v1 + p_v2) / 2
    JS = w1 * get_KL(p_v1, p_avg, 1) + w2 * get_KL(p_v2, p_avg, 1)
    return JS


def get_exp_isp(dict_list, max_int):  ## 提取实验谱中的所有同位素峰
    exp_isp = []
    for i in range(len(dict_list)):
        dict_isp = dict_list[i]
        for x in dict_isp.keys():
            tmp_mz = dict_isp[x][0]
            tmp_int = []
            for j in range(len(dict_isp[x][1])):
                tmp_int.append(dict_isp[x][1][j] / max_int * 100)
            exp_isp.append([tmp_mz, np.round(tmp_int, 4)])
    return exp_isp


def get_relative_sp(exp_sp, max_int):
    sp_with_rela_int = []
    for x in exp_sp:
        mass = x[0]
        rela_int = x[1] / max_int * 100
        sp_with_rela_int.append([mass, rela_int])
    return sp_with_rela_int


def sort_exp(exp_isp):
    exp_isp_copy = exp_isp.copy()
    for i in range(len(exp_isp_copy)):
        j = i
        min_i = i
        min_left = exp_isp_copy[i][0][0]
        while (j < len(exp_isp_copy)):
            if exp_isp_copy[j][0][0] < min_left:
                min_i = j
                min_left = exp_isp_copy[j][0][0]
            j = j + 1

        temp = exp_isp_copy[min_i]
        exp_isp_copy[min_i] = exp_isp_copy[i]
        exp_isp_copy[i] = temp
    return exp_isp_copy


def sort_the(the_isp, the_lost_list, the_z_list, the_H_Na):
    the_isp_copy = the_isp.copy()
    the_lost_list_copy = the_lost_list.copy()
    the_z_list_copy = the_z_list.copy()
    the_H_Na_copy = the_H_Na.copy()
    for i in range(len(the_isp_copy)):
        j = i + 1
        min_i = i
        min_left = the_isp_copy[i][0][0]
        while (j < len(the_isp_copy)):
            if the_isp_copy[j][0][0] < min_left:
                min_i = j
                min_left = the_isp_copy[j][0][0]
            j = j + 1

        temp_the = the_isp_copy[min_i]
        temp_lost = the_lost_list_copy[min_i]
        temp_z = the_z_list_copy[min_i]
        temp_hna = the_H_Na_copy[min_i]

        the_isp_copy[min_i] = the_isp_copy[i]
        the_lost_list_copy[min_i] = the_lost_list_copy[i]
        the_z_list_copy[min_i] = the_z_list_copy[i]
        the_H_Na_copy[min_i] = the_H_Na_copy[i]

        the_isp_copy[i] = temp_the
        the_lost_list_copy[i] = temp_lost
        the_z_list_copy[i] = temp_z
        the_H_Na_copy[i] = temp_hna
    return the_isp_copy, the_lost_list_copy, the_z_list_copy, the_H_Na_copy


def get_high_peak(all_mass_list, dict_match_exp):
    result = []
    max_int = []
    for x in all_mass_list:
        tmp_max = 0
        tmp_idx = 0
        for i in range(len(x)):
            mz = x[i]
            tmp_int = dict_match_exp[mz]
            if tmp_int > tmp_max:
                tmp_max = tmp_int
                tmp_idx = i
        result.append(x[tmp_idx])
        max_int.append(tmp_max)
    return result, max_int


def get_total_intensity(filter_mz, max_int):
    total_int = 0
    for peak in filter_mz:
        total_int += peak[1]
    total_int = np.round(total_int / max_int * 100, 4)
    return total_int

