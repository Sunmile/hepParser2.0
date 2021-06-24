import numpy as np
import pandas as pd
import copy
import pickle as pk
import os
import matplotlib.pyplot as plt
from src.enumerate import *
from src.utils import *
from brainpy import isotopic_variants


def get_one_Z_isotopic(MZ_list, fit_list, dict_list, ppm=10, Z=1, isotopic_num=5, js_sim_th=0.0):
    fit_list = fit_list[0:isotopic_num]
    max_int = 0
    for i in range(len(MZ_list)):
        if MZ_list[i][1] > max_int:
            max_int = MZ_list[i][1]
    visited_list = np.zeros(len(MZ_list))
    isotopic_record = np.zeros((len(MZ_list), isotopic_num))
    tmp_com = {}
    # if Z == 2:
    #     tmp_com = dict_list[0]
    # if Z == 4:
    #     tmp_com = dict_list[1]
    dict = {}
    dict_theo = {}
    dict_cos = {}
    for i in range(len(MZ_list)):
        if visited_list[i] == 0:
            M = MZ_list[i][0]
            flag_M = M
            flag_num = 1
            mass_l = [M]
            isotopic_record[i][0] = MZ_list[i][1]
            visited_list[i] = 1
            # theory_iso_list = generate_isotope_peak(M, Z, isotopic_num)
            j = i
            tmp_vis = []
            while j in range(i, len(MZ_list)):
                Mj = MZ_list[j][0]
                win = Mj * ppm * 0.000001
                dist = Mj - flag_M
                if dist < (1.0 / Z - win):
                    j += 1
                elif dist <= (1.0 / Z + win):
                    isotopic_record[i][flag_num] = MZ_list[j][1]
                    mass_l.append(Mj)
                    tmp_vis.append(j)
                    flag_M += 1.0 / Z
                    flag_num += 1
                    if flag_num >= isotopic_num:
                        break
                    j += 1
                else:
                    flag_M += 1.0 / Z
                    mass_l.append(mass_l[-1] + 1 / Z)
                    flag_num += 1
                    if flag_num >= isotopic_num:
                        break
            if np.count_nonzero(isotopic_record[i]) >= 3 \
                    and np.sum(isotopic_record[i]) >= max_int * 0.005 \
                    and isotopic_record[i][1] + isotopic_record[i][3] > 0:
                theory_distribution = [x(M * Z) if x(M * Z) > 0 else 0 for x in fit_list]
                tmp_js_sim = 1 - get_JS(theory_distribution, isotopic_record[i])
                if tmp_js_sim > js_sim_th:
                    if M in tmp_com.keys():
                        tmp_com.pop(M)
                    dict[M] = [mass_l, isotopic_record[i]]
                    dict_theo[M] = [mass_l, theory_distribution]
                    dict_cos[M] = tmp_js_sim
                    # for x in tmp_vis:
                    #     visited_list[x] = 1
    # if Z == 2:
    #     dict_list[0] = tmp_com
    # if Z == 4:
    #     dict_list[1] = tmp_com
    return dict_cos, dict, dict_theo


# 从过滤后的质谱中找所有的同位素峰，返回每个价态的同位素峰的字典，以及所有同位素峰的合并list
def get_isotopic(MZ_list, fit_list, ppm=10, max_Z=7, isotopic_num=4, js_sim_th=0.0):
    dict_list = []
    dict_score_list = []
    dict_the_iso_list = []
    for z in range(1, max_Z + 1):
        dict_score, dict_original, dict_the_iso = get_one_Z_isotopic(MZ_list, fit_list, dict_list, ppm, z, isotopic_num,
                                                                     js_sim_th)

        dict_list.append(dict_original)
        dict_score_list.append(dict_score)
        dict_the_iso_list.append(dict_the_iso)
        # print('Z='+str(z))
    return dict_score_list, dict_list, dict_the_iso_list


def get_all_isotope_distribute(filter_list):
    dict_all_mass = {}
    all_mass = []
    result = []
    for i in range(len(filter_list)):
        tmp_mass = get_molecular_mass(filter_list[i])
        all_mass.append(tmp_mass)
        dict_all_mass[tmp_mass] = filter_list[i]
    all_mass = sorted(all_mass)
    flag = all_mass[0]
    for x in all_mass:
        if 100 <= x <= 4000 and x - flag > 5:
            tmp_com = dict_all_mass[x]
            tmp_lwh = {}
            tmp_lwh['H'] = tmp_com[0]
            tmp_lwh['C'] = tmp_com[1]
            tmp_lwh['N'] = tmp_com[2]
            tmp_lwh['O'] = tmp_com[3]
            tmp_lwh['S'] = tmp_com[4]
            tmp_distribution = isotopic_variants(tmp_lwh, npeaks=5, charge=0)
            tmp_point = []
            tmp_point.append(x)
            for peak in tmp_distribution:
                tmp_point.append(peak.intensity)
            result.append(tmp_point)
            flag = x
    return result


# 拟合同位素中 峰的高度，得到intensity随mass的变化函数
def fit_all_point(point_list, degree=2):
    x = []
    y1 = []
    y2 = []
    y3 = []
    y4 = []
    y5 = []
    for i in point_list:
        x.append(i[0])
        y1.append(i[1])
        y2.append(i[2])
        y3.append(i[3])
        y4.append(i[4])
        y5.append(i[5])
    z1 = np.polyfit(x, y1, degree)
    z2 = np.polyfit(x, y2, degree)
    z3 = np.polyfit(x, y3, degree)
    z4 = np.polyfit(x, y4, degree)
    z5 = np.polyfit(x, y5, degree)
    p1 = np.poly1d(z1)
    p2 = np.poly1d(z2)
    p3 = np.poly1d(z3)
    p4 = np.poly1d(z4)
    p5 = np.poly1d(z5)
    yvals1 = p1(x)
    yvals2 = p2(x)
    yvals3 = p3(x)
    yvals4 = p4(x)
    yvals5 = p5(x)
    plt.figure(figsize=(10, 8))
    plt.scatter(x, y1, c='r')
    plt.plot(x, yvals1, 'r')
    plt.scatter(x, y2, c='g')
    plt.plot(x, yvals2, 'g')
    plt.scatter(x, y3, c='b')
    plt.plot(x, yvals3, 'b')
    plt.scatter(x, y4, c='k')
    plt.plot(x, yvals4, 'k')
    plt.scatter(x, y5, c='y')
    plt.plot(x, yvals5, 'y')
    plt.savefig('fit.png')
    plt.show()
    return [p1, p2, p3, p4, p5]


def get_fit_pk(dir, dp):
    dir = dir + str(dp) + '.pk'
    if os.path.exists(dir):
        with open(dir, 'rb') as f:
            data = pk.load(f)
        return data
    else:
        # all_list = enumerate_component(20)  # 枚举所有可能肝素分子的构成
        # dict_dp_list = enumerate_dalteparin(dp)
        dict_dp_list = enumerate_enoxaparin(dp)
        all_list = dict_dp_list[dp]
        all_list = filter_component(all_list, 100, 5000)
        dict_mass_comp, dict_mass_atom, all_comp_list, all_atoms_list, all_mass_list = transform_component_to_atom(
            all_list)
        point_list = get_all_isotope_distribute(all_atoms_list)
        fit_list = fit_all_point(point_list, 4)
        save_pk(fit_list, dir)
        return fit_list


def get_fit_pk(dir):
    dir = dir + '.pk'
    if os.path.exists(dir):
        with open(dir, 'rb') as f:
            data = pk.load(f)
        return data
    else:
        all_list = enumerate_component(20)  # 枚举所有可能肝素分子的构成
        # dict_dp_list = enumerate_dalteparin(dp)
        # dict_dp_list = enumerate_enoxaparin(dp)
        # all_list = dict_dp_list[dp]
        all_list = filter_component(all_list, 100, 5000)
        dict_mass_comp, dict_mass_atom, all_comp_list, all_atoms_list, all_mass_list = transform_component_to_atom(all_list)
        point_list = get_all_isotope_distribute(all_atoms_list)
        fit_list = fit_all_point(point_list, 4)
        save_pk(fit_list, dir)
        return fit_list