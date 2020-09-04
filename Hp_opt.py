import numpy as np
import pandas as pd
import copy
import itertools
from pprint import pprint as pl
import pickle as pk
import os
from time import time
import math
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from sklearn.linear_model import Lasso
from brainpy import isotopic_variants
from collections import Counter as ct
from peak_align import *

dict_atom = {}
dict_atom['H1'] = 1.00783
dict_atom['H2'] = 2.01410
dict_atom['C12'] = 12.00000
dict_atom['C13'] = 13.00335
dict_atom['N14'] = 14.00307
dict_atom['N15'] = 15.00011
dict_atom['O16'] = 15.99492
dict_atom['O17'] = 16.99913
dict_atom['O18'] = 17.99916
dict_atom['S32'] = 31.97207
dict_atom['S33'] = 32.97146
dict_atom['S34'] = 33.96787
dict_atom['S36'] = 35.96708
dict_atom['Na'] = 22.98977


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
            if MZ >= 350:
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
        if MZ >= 350:
            if intensity > max_int:
                max_int = intensity
            MZ_int.append([MZ, intensity])
    for i in range(len(MZ_int)):
        reletive_int = MZ_int[i][1] / max_int
        MZ_int[i].append(np.round(reletive_int, 6))
    return MZ_int, max_int


def filter_spectra(MZ_int, th, max_int, min_int=500):
    filter_MZ = []
    for i in range(len(MZ_int)):
        if MZ_int[i][1] > min_int and MZ_int[i][2] > th:
            filter_MZ.append(MZ_int[i])
    return filter_MZ


def save_file(data, dir):
    f = open(dir, 'w', newline='')
    for i in range(len(data)):
        f.writelines(str(data[i][0]) + ' ' + str(data[i][1]))
        f.writelines('\n')
    f.close()


def generate_isotope_peak(M, Z=1, isotopic_num=4):
    iso_list = []
    for j in range(isotopic_num):
        iso_list.append(M + j / float(Z))
    return iso_list


def get_distance_list(MZ_list):
    distance_list = []
    for j in range(1, len(MZ_list)):
        distance_list.append(MZ_list[j][0] - MZ_list[j - 1][0])
    return distance_list


def get_one_Z_isotopic(MZ_list, fit_list, dict_list, ppm=10, Z=1, isotopic_num=5, cos_sim_th=0):
    fit_list = fit_list[0:isotopic_num]
    visited_list = np.zeros(len(MZ_list))
    isotopic_record = np.zeros((len(MZ_list), isotopic_num))
    tmp_com = {}
    if Z == 2:
        tmp_com = dict_list[0]
    if Z == 4:
        tmp_com = dict_list[1]
    dict = {}
    dict_filter = {}
    dict_cos = {}
    # if Z ==1 or Z==2:
    #     ppm = ppm*2
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
                    and np.sum(isotopic_record[i]) >= 10000 \
                    and isotopic_record[i][1] + isotopic_record[i][3] > 0:
                theory_distribution = [x(M * Z) if x(M * Z) > 0 else 0 for x in fit_list]
                tmp_cos_sim = 1 - get_JS(theory_distribution, isotopic_record[i])
                if tmp_cos_sim > cos_sim_th:
                    if M in tmp_com.keys():
                        tmp_com.pop(M)
                    dict[M] = [mass_l, isotopic_record[i]]
                    for x in tmp_vis:
                        visited_list[x] = 1
    if Z == 2:
        dict_list[0] = tmp_com
    if Z == 4:
        dict_list[1] = tmp_com
    # for m in dict.keys():
    #     test_mass_list, test_Int_list = dict[m]
    #     max_cos_sim = 0
    #     max_cos_pos = -2
    #     max_cos_Int_list = []
    #     max_cos_theory_list = []
    #     for start_p in range(-2, 3):
    #         tmp_m = m + start_p
    #         theory_distribution = [x(tmp_m * Z) if x(tmp_m * Z) > 0 else 0 for x in fit_list]
    #         test_distribution = [test_Int_list[i] if 0 <= i < isotopic_num
    #                              else 0 for i in range(start_p, isotopic_num + start_p)]
    #         tmp_cos_sim = cos_sim(theory_distribution, test_distribution)
    #         if tmp_cos_sim > max_cos_sim:
    #             max_cos_sim = tmp_cos_sim
    #             max_cos_pos = start_p
    #             max_cos_Int_list = test_distribution.copy()
    #             max_cos_theory_list = theory_distribution.copy()
    #     if max_cos_sim > cos_sim_th:
    #         tmp_m = m + max_cos_pos
    #         dict_filter[tmp_m] = max_cos_Int_list
    #         dict_cos[tmp_m] = np.round(max_cos_sim, 5)
    # print(test_mass_list, np.round(max_cos_sim, 5), max_cos_Int_list, np.around(max_cos_theory_list, decimals=4))
    return dict_filter, dict_cos, dict


# 从过滤后的质谱中找所有的同位素峰，返回每个价态的同位素峰的字典，以及所有同位素峰的合并list
def get_isotopic(MZ_list, fit_list, ppm=10, max_Z=5, isotopic_num=4, cos_sim_th=0):
    dict_list = []
    dict_change_start_list = []
    dict_cos_list = []
    isotopic_list_z = []
    for z in range(1, max_Z + 1):
        isotopic_dict, dict_cos, dict_original = get_one_Z_isotopic(MZ_list, fit_list, dict_list, ppm, z, isotopic_num,
                                                                    cos_sim_th)
        dict_change_start_list.append(isotopic_dict)
        dict_list.append(dict_original)
        dict_cos_list.append(dict_cos)
        # print('Z='+str(z))
        iso_list = []
        for x in isotopic_dict.keys():
            if x not in iso_list:
                iso_list.append(x)
            # print(str(x)+': '+ str(isotopic_dict[x]))
        isotopic_list_z.append(iso_list)
    return dict_change_start_list, dict_cos_list, dict_list


# 根据分子式计算分子质量
def get_molecular_mass(atom_list):
    Mass = 0
    atom_mass = [dict_atom['H1'], dict_atom['C12'], dict_atom['N14'],
                 dict_atom['O16'], dict_atom['S32'], dict_atom['Na']]
    for i in range(len(atom_list)):
        Mass += atom_mass[i] * atom_list[i]
    return np.round(Mass, 5)


# 考虑所有可能的肝素组成，返回所有肝素组成式

def enumerate_component(n=20):
    all_list = []  # bubaohe, baohe, putaotangan, yixian, liusuan, neimi, ganlu
    for i in range(0, n):  # i 糖醛酸数目
        if i - 1 <= 0:
            num_bao = [i, i + 1]
        else:
            num_bao = [i - 1, i, i + 1]
        for j in num_bao:  # j 葡萄糖胺数目
            if j == i + 1:
                num_bu = [0, 1]
                num_nei = [0]
                num_gan = [0, 1]
            elif j == i:
                num_bu = [0, 1]
                num_nei = [0, 1]
                num_gan = [0]
            else:
                num_bu = [0]
                num_nei = [0, 1]
                num_gan = [0]
            for x in num_bu:
                for y in num_nei:
                    for z in num_gan:
                        for m in range(j):
                            for n in range(3 * j - m + i + x + y + z):
                                tmp_comp = [x, i, j, m, n, y, z]
                                all_list.append(tmp_comp)
    print('All components counts', len(all_list))
    return all_list


def get_component_mass(one_comp):
    comp_atoms = [
        [8, 6, 0, 6, 0],  # 不饱和糖醛酸 HexA
        [10, 6, 0, 7, 0],  # 饱和糖醛酸  GlcA
        [13, 6, 1, 5, 0],  # 葡萄糖胺   GlcN
        [4, 2, 0, 2, 0],  # 乙酰基    Ac
        [2, 0, 0, 4, 1],  # 硫酸基    SO3
        [11, 6, 1, 4, 0],  # 内醚糖    Levoglucosan
        [12, 6, 0, 5, 0]  # 甘露糖    Man
    ]
    comp_atoms = np.array(comp_atoms)
    tmp_atoms = np.array([2, 0, 0, 1, 0])  # 带一个水
    H2O = np.array([2, 0, 0, 1, 0])  # 任两个基团结合，都脱一个水
    for i in range(len(one_comp)):
        tmp_atoms = tmp_atoms + comp_atoms[i] * one_comp[i] - H2O * one_comp[i]
    mass = get_molecular_mass(tmp_atoms)
    tmp_atoms = tmp_atoms.tolist()
    return np.round(mass, 5), tmp_atoms


def filter_component(all_list, lowth, highth):
    result_list = []
    for x in all_list:
        mass, tmp = get_component_mass(x)
        if lowth <= mass <= highth:
            result_list.append(x)
    return result_list


def transform_component_to_atom(all_list):
    dict_mass_comp = {}
    dict_mass_atom = {}
    all_atoms_list = []
    all_mass_list = []
    count = 0
    for x in all_list:
        mass, tmp_atoms = get_component_mass(x)
        all_atoms_list.append(tmp_atoms)
        all_mass_list.append(mass)
        if mass in dict_mass_atom.keys():
            print('There are two candidates with same mass!')
            print(x)
            print(dict_mass_comp[mass])
            print(tmp_atoms)
            print(dict_mass_atom[mass])

            dict_mass_atom[mass].append(tmp_atoms)
            dict_mass_comp[mass].append(x)

        else:
            dict_mass_comp[mass] = [x]
            dict_mass_atom[mass] = [tmp_atoms]
        count += 1
        print(count, mass, x)
    return dict_mass_comp, dict_mass_atom, all_list, all_atoms_list, all_mass_list


def generate_new_composition(filter_list, max_lost_count=20):
    dict_com = {}
    lost_record = []  # SO3, NH, NHSO3, COO, 丢失数目
    lost_list = [[0, 0, 0, 0, 0],  # 不丢失
                 [0, 0, 0, 3, 1],  # 丢失SO3
                 [1, 0, 1, 0, 0],  # 丢失NH
                 [1, 0, 1, 3, 1],  # 丢失NHSO3
                 [0, 1, 0, 2, 0]]  # 丢失COO
    all_lost_list = []
    ll = len(filter_list)
    for i in range(len(filter_list)):
        one_com = filter_list[i]
        tmp_list = []
        it_5_count = itertools.combinations_with_replacement([0, 1, 2, 3, 4], max_lost_count)
        for one_iter in it_5_count:
            count_dict = ct(one_iter)
            tmp_lost = [lost_list[x] for x in one_iter]
            sum_lost = np.sum(tmp_lost, axis=0).tolist()
            tmp_com = list(map(lambda x: x[0] - x[1], zip(one_com, sum_lost)))
            flag_list = [x >= 0 for x in tmp_com]
            if False in flag_list:
                continue
            # elif tmp_com not in all_lost_list:
            else:
                tmp_tup = tuple(tmp_com)
                if tmp_tup in dict_com.keys():
                    continue
                else:
                    dict_com[tmp_tup] = 0
                    lost_record.append([count_dict[1], count_dict[2], count_dict[3], count_dict[4]])
                    tmp_list.append(tmp_com)
        all_lost_list.extend(tmp_list)
    #     print(i, ll)
    # print(len(all_lost_list))
    return all_lost_list, lost_record


# 返回两个list，用于标注同位素峰，以及同位素最高峰价态
def get_isotope_results(dict_peak_list):
    result1 = []
    result2 = []
    tmp_csv = []
    dict_z = {}
    dict_int = {}
    dict_pos = {}
    for z in range(len(dict_peak_list)):
        dict_peaks = dict_peak_list[z]
        for mz in dict_peaks.keys():
            int_list = dict_peaks[mz]
            for j in range(len(int_list)):
                result1.append([np.round(mz + j / float(z + 1), 4), int_list[j]])
            pos = np.argmax(int_list)
            tmp_mass = mz + pos / float(z + 1)
            tmp_int = int_list[pos]
            tmp_csv.append([np.round(tmp_mass, 4), tmp_int, str(z + 1) + '-', pos + 1])
            if tmp_mass in dict_z.keys():
                dict_z[tmp_mass] += ',' + str(z + 1) + '-'
                dict_pos[tmp_mass] += ',' + str(pos + 1)
            else:
                dict_z[tmp_mass] = str(z + 1) + '-'
                dict_pos[tmp_mass] = str(pos + 1)
                dict_int[tmp_mass] = tmp_int
    for x in sorted(dict_z.keys()):
        tmp_int = dict_int[x]
        tmp_z = dict_z[x]
        tmp_pos = dict_pos[x]
        result2.append([np.round(x, 4), tmp_int, tmp_z, tmp_pos])
    tmp_pd = pd.DataFrame(tmp_csv)
    tmp_pd.to_csv('Isotopic_peak_position.csv', index=None, header=None)
    return result1, result2


def save_pk(data, dir):
    with open(dir, 'wb') as f:
        pk.dump(data, f)


def get_comp_pk(dir):
    if os.path.exists(dir):
        with open(dir, 'rb') as f:
            data = pk.load(f)
        return data
    else:
        all_list = enumerate_component(20)  # 枚举所有可能肝素分子的构成
        all_list = filter_component(all_list, 350, 3600)
        dict_mass_comp, dict_mass_atom, all_comp_list, all_atoms_list, all_mass_list = transform_component_to_atom(
            all_list)
        result = [dict_mass_comp, dict_mass_atom, all_comp_list, all_atoms_list, all_mass_list]
        save_pk(result, dir)
        return result


def get_fit_pk(dir):
    if os.path.exists(dir):
        with open(dir, 'rb') as f:
            data = pk.load(f)
        return data
    else:
        all_list = enumerate_component(20)  # 枚举所有可能肝素分子的构成
        all_list = filter_component(all_list, 350, 3600)
        dict_mass_comp, dict_mass_atom, all_atoms_list, all_mass_list = transform_component_to_atom(all_list)
        point_list = get_all_isotope_distribute(all_atoms_list)
        fit_list = fit_all_point(point_list, 4)
        save_pk(fit_list, dir)
        return fit_list


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
        if 350 <= x <= 3600 and x - flag > 5:
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


def generate_theory_SP(compsition, prob=0.5):
    all_com, lost_record = generate_new_composition([compsition], 2)
    # pl(all_com)
    lost_comp = []
    Z_list = []
    all_mz = []
    # charges = [-1]
    charges = [-1, -2, -3, -4, -5]
    for i in range(len(all_com)):
        x = all_com[i]
        losted = lost_record[i]
        tmp_com = x
        tmp_lwh = {}
        tmp_lwh['H'] = tmp_com[0]
        tmp_lwh['C'] = tmp_com[1]
        tmp_lwh['N'] = tmp_com[2]
        tmp_lwh['O'] = tmp_com[3]
        tmp_lwh['S'] = tmp_com[4]
        for j in charges:
            tmp_distribution = isotopic_variants(tmp_lwh, npeaks=5, charge=j)
            tmp_mz = []
            tmp_int = []
            for peak in tmp_distribution:
                tmp_mz.append(np.round(peak.mz, 4))
                tmp_int.append(peak.intensity * 100 * np.power(prob, sum(losted)))
                # all_mz.append([np.round(peak.mz,4), peak.intensity])
            all_mz.append([tmp_mz, np.round(tmp_int, 4)])
            lost_comp.append(losted)
            Z_list.append(-j)
    # sorted_mz = sorted(all_mz, key=lambda mz: mz[0])
    return all_mz, lost_comp, Z_list


def get_exp_isp(dict_list, max_int):  ## 提取实验谱中的所有同位素峰
    exp_isp = []
    for i in range(len(dict_list)):
        dict_isp = dict_list[i]
        for x in dict_isp.keys():
            tmp_mz = dict_isp[x][0]
            tmp_int = []
            for j in range(5):
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


def get_one_mass_comp(mass, z, ppm, the_HP):
    atom_list = the_HP[1]
    comp_mass = mass * z + z
    win = comp_mass * ppm * 0.000001
    result = []
    for m in atom_list.keys():
        one_atom = atom_list[m][0]
        all_list, lost_record = generate_new_composition([one_atom], 5)
        for i in range(len(all_list)):
            x = all_list[i]
            losted = lost_record[i]
            mass = get_molecular_mass(x)
            if abs(mass - comp_mass) <= win:
                comp = the_HP[0][m][0]
                result.append([comp, losted])
                print(comp_mass, mass)
                print(comp, losted)
    return result


def get_JS(v1, v2, w1=0.5, w2=0.5):
    v1 = np.array(v1)
    v2 = np.array(v2)

    p_v1 = v1 / sum(v1)
    p_v2 = v2 / sum(v2)
    p_avg = (p_v1 + p_v2) / 2
    JS = w1 * get_KL(p_v1, p_avg, 1) + w2 * get_KL(p_v2, p_avg, 1)
    return JS


def cal_all_the_sp(all_atom_list, prob):
    all_the = []
    all_lost = []
    all_z = []
    all_the_01 = []
    for one_atom in all_atom_list:
        the_isp, lost_comp, Z_list = generate_theory_SP(one_atom, prob)
        the_isp, lost_comp, Z_list = sort_the(the_isp, lost_comp, Z_list)
        min_id = 0
        for i in range(len(the_isp)):
            x = the_isp[i][0][0]
            if x >= 350:
                min_id = i
                break
        the_isp = the_isp[min_id:]
        lost_comp = lost_comp[min_id:]
        Z_list = Z_list[min_id:]
        all_the.append(np.float32(the_isp))
        # all_the.append(the_isp)
        all_lost.append(lost_comp)
        all_z.append(Z_list)
        one_01 = change_the_format(the_isp)
        all_the_01.append(one_01)
    result = [all_the, all_lost, all_z, all_the_01]
    return result


def get_the_sp(dir, all_atom_list, prob):
    if os.path.exists(dir):
        with open(dir, 'rb') as f:
            data = pk.load(f)
        return data
    else:
        data = cal_all_the_sp(all_atom_list, prob)
        save_pk(data, dir)
        return data


def change_the_format(sp, win=0.1):
    index_t = []
    min_th = 350 * np.int(1 / win)
    max_th = 1200 * np.int(1 / win)
    for one_isp in sp:
        for i in range(len(one_isp[0])):
            mz = one_isp[0][i]
            int = one_isp[1][i]
            if int > 0:
                mass = np.int(mz * np.int(1 / win))
                if min_th <= mass < max_th:
                    index_t.append(mass - min_th)
    return index_t


def transform_the_01(sp, win=0.1):
    result = np.zeros([850 * np.int(1 / win)], dtype=bool)
    result[sp] = True
    return result


def change_sp_format(sp, win=0.05):
    result = np.zeros([850 * np.int(1 / win)], dtype=bool)
    min_th = 350 * np.int(1 / win)
    max_th = 1200 * np.int(1 / win)
    for one_isp in sp:
        for i in range(len(one_isp[0])):
            mz = one_isp[0][i]
            int = one_isp[1][i]
            if int > 0:
                mass = np.int(mz * np.int(1 / win))
                if min_th <= mass < max_th:
                    result[mass - min_th] = True
    return result


def compared_score(a, b):
    # return len([1 for i, j in zip(list(a), list(b)) if i==j==True])
    # x = sum((b==True)&(a==True))
    x = np.count_nonzero((b == True) & (a == True))
    return x


def score_match_isp(isp_1, isp_2, ppm):
    '''
    input:
    isp_1: 同位素峰簇1   理论
    isp_2: 同位素峰簇2   实际
    example: [[m/z1,m/z2,m/z3,m/z4,m/z5],[int1,int2,int3,int4,int5]]
    注意: isp_1, isp_2 的m/z 值必须对应相等
    output:

    score: 匹配得分
    '''
    mz_list1 = isp_1[0]
    mz_list2 = isp_2[0]
    int_list1 = isp_1[1]
    int_list2 = isp_2[1]
    int_weight = 0
    for x in int_list2:
        int_weight += np.log(x + 1)
        # int_weight += x
    # int_weight = np.log(np.sum(int_list2)+1)
    score_mz = 0
    for i in range(len(mz_list1)):
        diff = abs(mz_list1[i] - mz_list2[i])
        win = ppm * mz_list1[i] * 0.000001
        if diff <= win:
            score_mz += diff / win + 0.1
    # score_int = cos_sim(int_list1,int_list2)
    score_int = int_weight * (1 - get_JS(int_list1, int_list2, 0.5, 0.5))
    # score_int = get_KL(int_list1,int_list2)
    # score_all = score_mz * score_int
    score_all = score_int
    return score_all


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


def sort_the(the_isp, the_lost_list, the_z_list):
    the_isp_copy = the_isp.copy()
    the_lost_list_copy = the_lost_list.copy()
    the_z_list_copy = the_z_list.copy()
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

        the_isp_copy[min_i] = the_isp_copy[i]
        the_lost_list_copy[min_i] = the_lost_list_copy[i]
        the_z_list_copy[min_i] = the_z_list_copy[i]

        the_isp_copy[i] = temp_the
        the_lost_list_copy[i] = temp_lost
        the_z_list_copy[i] = temp_z
    return the_isp_copy, the_lost_list_copy, the_z_list_copy


def match_isotopic(exp_isp, the_isp, ppm, thresh, min_match_num, the_lost_list, the_z_list, p=0.9):
    '''
    input:
    exp_isp: 实验谱，形式是同位素峰簇的集合，每个同位素峰簇包含五个峰的m/z值，以及对应的int值。三重list，
    example: [[[m/z1,m/z2,m/z3,m/z4,m/z5],[int1,int2,int3,int4,int5]],[[m/z1,m/z2,m/z3,m/z4,m/z5],[int1,int2,int3,int4,int5]]]
    the_isp: 理论谱，格式与exp_isp 相同
    ppm: win = mass*ppm*10^-6
    win: 匹配窗口的大小，两个峰m/z 差值在win之内，被认为是匹配成功
    thresh: 单个同位素峰匹配得分阈值
    min_match_num: 	单个同位素峰簇中至少匹配的峰数


    output:
    global_score: 匹配的全局得分, 单次匹配的加和
    matched_isp_list: 被匹配上的同位素峰簇的集合
    ret_the_isp_list: 匹配上的同位素峰簇的对应理论峰簇的集合
    score_list: 单个同位素峰簇得分list
    score_mass: 匹配的单同位素峰

    '''

    '''
    执行操作:
    对于exp_isp 中的每一个 同位素峰簇，分别与 the_isp 中的同位素峰簇匹配
    若有不低于min_match_num个峰被匹配上，则认为该两个同位素峰簇匹配上了。
    匹配上的话，以 the_isp中的同位素峰簇为准，将exp_isp中的同位素峰簇补全，补充的峰的int值为0(注意，这步不能修改原exp_isp值),
    然后调用score_match_isp(isp_1,isp2,win) 获得 匹配得分。
    example1：
    exp_isp中：isp_1=[[1, 1.5, 2, 2.5, 3],[10, 20, 11, 10, 5]]
    the_isp中：isp_2=[[1, 2, 3, 4, 5],[20, 5, 10, 22, 4]]
    则 isp_1 更改为：[[1, 2, 3, 4, 5],[10, 11, 5, 0, 0]]
    example2：
    exp_isp中：isp_1=[[3, 4, 5, 6, 7],[10, 20, 11, 10, 5]]
    the_isp中：isp_2=[[1, 2, 3, 4, 5],[20, 5, 10, 22, 4]]
    则 isp_1 更改为：[[1, 2, 3, 4, 5],[0, 0, 10, 20, 11]]

    只有匹配得分超过thresh的，才会被加到matched_isp_list中, 加的是更改后的exp_isp 中的同位素峰簇。
    '''
    # 找到最长的间距
    MAX_L = 0
    for item in exp_isp:
        temp = item[0][-1] - item[0][0]
        if temp > MAX_L:
            MAX_L = temp

    for item in the_isp:
        temp = item[0][-1] - item[0][0]
        if temp > MAX_L:
            MAX_L = temp

    matched_isp_list = []
    matched_score_list = []
    global_score = 0
    # 记录matched_isp_list中每组匹配的实验谱和理论普索引，用于后期筛选
    index_exp = []
    index_the = []
    i = 0
    j0 = 0
    while (i < len(exp_isp)):
        # 对于每个实验谱中的峰簇，找到匹配分数最大的理论谱
        exp_item = exp_isp[i]
        while (exp_item[0][0] > the_isp[j0][0][0] + MAX_L):
            j0 = j0 + 1
            if j0 == len(the_isp):
                j0 = j0 - 1
                break

        if j0 == len(the_isp) - 1 and exp_item[0][0] > the_isp[j0][0][0] + MAX_L:
            break

        max_score = 0
        max_list = None
        max_j = j0
        j = j0
        while (j < len(the_isp)):
            the_item = the_isp[j]
            tmp_lost_sum = np.sum(the_lost_list[j])
            if the_item[0][0] > exp_item[0][0] + MAX_L:
                break

            counter = 0
            # 计算两个峰簇的匹配数量,并整理出补全后的峰簇temp_exp
            temp_exp = []
            temp_exp.append([])
            temp_exp.append([])
            match_str = ''
            win = ppm * the_item[0][0] * 0.000001
            for aim_mz in the_item[0]:
                # 两个峰簇匹配时，使用理论谱峰簇的m/z值在实验谱中搜索，如果找到win之内的，则记录，否则记0
                temp_exp[0].append(aim_mz)
                temp_value = 0
                k = 0
                match_flag = False
                for k in range(counter, len(exp_item[0])):
                    exp_mz = exp_item[0][k]
                    if abs(exp_mz - aim_mz) <= win and exp_item[1][k] != 0:
                        temp_value = exp_item[1][k]
                        temp_exp[0][-1] = exp_item[0][k]
                        counter = counter + 1
                        match_str += '1'
                        match_flag = True
                        break
                temp_exp[1].append(temp_value)
                if not match_flag:
                    match_str += '0'
            # 匹配数目足够，则将该匹配看做备选项
            if counter < min_match_num or match_str == '10101':
                j = j + 1
                continue

            # 计算备选项匹配得分
            this_score = score_match_isp(the_item, temp_exp, ppm)
            if tmp_lost_sum > 0:
                this_score = p * this_score

            if this_score <= thresh:
                j = j + 1
                continue

            # 备选项的分数比当前最大值大的话，更新最大值及相关索引
            if this_score > max_score:
                max_score = this_score
                max_list = temp_exp.copy()
                max_j = j
            j = j + 1

        # 若该实验峰簇未找到匹配项
        if max_score == 0:
            i = i + 1
            continue

        # 记录补全后的峰簇、累加匹配分数、记录两个索引
        matched_isp_list.append(max_list)
        matched_score_list.append(max_score)
        global_score = global_score + max_score
        index_exp.append(i)
        index_the.append(max_j)
        i = i + 1

    # 判断是否有重复
    if len(set(index_the)) == len(index_the):
        ret_the_isp_list = []
        ret_the_lost_list = []
        ret_the_z_list = []
        for i in index_the:
            ret_the_isp_list.append(the_isp[i])
            ret_the_lost_list.append(the_lost_list[i])
            ret_the_z_list.append(the_z_list[i])
        score_mass = []  # 标记 得分的单同位素峰
        for i in range(len(matched_isp_list)):
            mass = matched_isp_list[i][0]
            score_mass.append(mass)
        return global_score, matched_isp_list, ret_the_isp_list, matched_score_list, score_mass, ret_the_lost_list, ret_the_z_list

    # 有重复的话
    index = []
    # 找到重复的项，存入做为set的index中
    for i in range(len(index_the)):
        j = i + 1
        while (j < len(index_the)):
            if (index_the[i] == index_the[j]):
                index_flag = 0
                for item in index:
                    if i in item:
                        item.add(j)
                        index_flag = 1
                        break
                    elif j in item:
                        item.add(i)
                        index_flag = 1
                        break
                if index_flag == 0:
                    temp_set = set()
                    temp_set.add(i)
                    temp_set.add(j)
                    index.append(temp_set)
            j = j + 1

    delete_counter = 0
    # 对于每一组重复
    for item in index:
        item_list = list(item)
        if len(item) == 2:
            i = item_list[0]
            j = item_list[1]
            # flag表征两类情况，1代表可合并，0代表不可合并应取分数最大
            flag = 0
            # if和elif代表两种可能出现合并的起始情况，若可合并，则最后flag会被置为1
            # if：i对应峰簇前半部分为0，j对应峰簇后半部分为0
            if matched_isp_list[i][1][0] == 0 and matched_isp_list[j][1][0] != 0:
                # 核心思路是找到i由0变非0，以及j由非0变0的位置，使用change做标记
                # 若上述改变只有一次，且未发生其他一些特殊情况，则认为i和j能够合并
                # elif里思路相同
                change = 0
                for k in range(len(matched_isp_list[i][1])):
                    if matched_isp_list[i][1][k] != 0 and matched_isp_list[j][1][k] != 0:
                        break
                    if change == 0 and matched_isp_list[i][1][k] == 0 and matched_isp_list[j][1][k] != 0:
                        continue
                    elif change == 0 and matched_isp_list[j][1][k] == 0 and matched_isp_list[i][1][k] == 0:
                        if k == len(matched_isp_list[i][1]) - 1:
                            flag = 1
                        continue
                    elif change == 0 and matched_isp_list[j][1][k] == 0 and matched_isp_list[i][1][k] != 0:
                        change = 1
                        if k == len(matched_isp_list[i][1]) - 1:
                            flag = 1
                        continue
                    elif change == 1 and matched_isp_list[i][1][k] == 0 and matched_isp_list[j][1][k] == 0:
                        if k == len(matched_isp_list[i][1]) - 1:
                            flag = 1
                        continue
                    elif change == 1 and matched_isp_list[i][1][k] != 0 and matched_isp_list[j][1][k] == 0:
                        if k == len(matched_isp_list[i][1]) - 1:
                            flag = 1
                        continue
                    elif change == 1 and matched_isp_list[i][1][k] == 0 and matched_isp_list[j][1][k] != 0:
                        break

            # elif：i对应峰簇后半部分为0，j对应峰簇前半部分为0
            elif matched_isp_list[i][1][0] != 0 and matched_isp_list[j][1][0] == 0:
                change = 0
                for k in range(len(matched_isp_list[i][1])):
                    if matched_isp_list[i][1][k] != 0 and matched_isp_list[j][1][k] != 0:
                        break
                    if change == 0 and matched_isp_list[i][1][k] != 0 and matched_isp_list[j][1][k] == 0:
                        continue
                    elif change == 0 and matched_isp_list[i][1][k] == 0 and matched_isp_list[j][1][k] == 0:
                        if k == len(matched_isp_list[i][1]) - 1:
                            flag = 1
                        continue
                    elif change == 0 and matched_isp_list[i][1][k] == 0 and matched_isp_list[j][1][k] != 0:
                        change = 1
                        if k == len(matched_isp_list[i][1]) - 1:
                            flag = 1
                    elif change == 1 and matched_isp_list[i][1][k] == 0 and matched_isp_list[j][1][k] == 0:
                        if k == len(matched_isp_list[i][1]) - 1:
                            flag = 1
                        continue
                    elif change == 1 and matched_isp_list[i][1][k] == 0 and matched_isp_list[j][1][k] != 0:
                        if k == len(matched_isp_list[i][1]) - 1:
                            flag = 1
                        continue
                    elif change == 1 and matched_isp_list[i][1][k] != 0 and matched_isp_list[j][1][k] == 0:
                        break

            # 如果可以合并
            if flag == 1:
                # 以i为模板制作一个新的峰簇，用j将其补全
                temp_isp = []
                temp_isp.append(matched_isp_list[i][0].copy())
                temp_isp.append(matched_isp_list[i][1].copy())
                for k in range(len(temp_isp[0])):
                    if matched_isp_list[i][1][k] != 0:
                        continue
                    else:
                        # 同时更新m/z信息和峰强信息
                        temp_isp[0][k] = matched_isp_list[j][0][k]
                        temp_isp[1][k] = matched_isp_list[j][1][k]

                # 计算分数，减去原来的，加上新的
                score1 = matched_score_list[i]
                score2 = matched_score_list[j]
                score3 = score_match_isp(the_isp[index_the[i]], temp_isp, ppm)
                if np.sum(the_lost_list[index_the[i]]) > 0:
                    score3 = p * score3
                global_score = global_score - score1 - score2 + score3
                # 将原来的峰簇标记
                matched_isp_list[i] = 'deleted'
                matched_isp_list[j] = 'deleted'
                matched_score_list[i] = -1
                matched_score_list[j] = -1
                index_the.append(index_the[i])
                index_the[i] = -1
                index_the[j] = -1
                # 记录应删除峰簇的数量
                delete_counter = delete_counter + 2
                # 加入合并后的峰簇
                matched_isp_list.append(temp_isp)
                matched_score_list.append(score3)

            # 无须合并，找匹配分数最大的进行保留
            else:
                score1 = matched_score_list[i]
                score2 = matched_score_list[j]
                if score1 < score2:
                    global_score = global_score - score1
                    matched_isp_list[i] = 'deleted'
                    matched_score_list[i] = -1
                    index_the[i] = -1
                else:
                    global_score = global_score - score2
                    matched_isp_list[j] = 'deleted'
                    matched_score_list[j] = -1
                    index_the[j] = -1
                delete_counter = delete_counter + 1

        else:
            temp_score_list = []
            temp_max_score = 0
            temp_index = -1
            i = 0
            for index_item in item_list:
                score = matched_score_list[index_item]
                temp_score_list.append(score)
                if score > temp_max_score:
                    temp_max_score = score
                    temp_index = i
                i = i + 1

            for j, score in enumerate(temp_score_list):
                if j == temp_index:
                    continue
                global_score = global_score - score
                matched_isp_list[item_list[j]] = 'deleted'
                matched_score_list[item_list[j]] = -1
                index_the[item_list[j]] = -1
                delete_counter = delete_counter + 1

    # 根据标记数量，删除标记峰簇
    for i in range(delete_counter):
        matched_isp_list.remove('deleted')
        index_the.remove(-1)
        matched_score_list.remove(-1)

    # 更新理论谱list
    ret_the_isp_list = []
    ret_the_lost_list = []
    ret_the_z_list = []
    for i in index_the:
        ret_the_isp_list.append(the_isp[i])
        ret_the_lost_list.append(the_lost_list[i])
        ret_the_z_list.append(the_z_list[i])

    score_mass = []  # 标记 得分的单同位素峰
    for i in range(len(matched_isp_list)):
        mass = matched_isp_list[i][0]
        score_mass.append(mass)

    return global_score, matched_isp_list, ret_the_isp_list, matched_score_list, score_mass, ret_the_lost_list, ret_the_z_list


def calculate_all_hp_score(the_spectra, dict_list, the_HP, max_int, ppm=10, thresh=0, min_match_num=3, prob=0.5):
    exp_isp = get_exp_isp(dict_list, max_int)
    exp_isp = sort_exp(exp_isp)
    exp_isp_01 = change_sp_format(exp_isp, 0.1)
    all_comp = []
    all_atom = []
    all_comp_mass = []
    all_comp_score = []
    all_comp_match = []
    all_comp_lost = []
    all_comp_z = []
    dict_match_exp = {}
    dict_theo_list = []
    all_score_list = []
    all_score_mass = []
    count = 0
    the_comp_list = the_HP[2]
    the_atom_list = the_HP[3]
    the_mass_list = the_HP[4]
    pass_count = 0
    ttt = time()
    for i in range(len(the_comp_list)):
        count += 1
        print("ID\t" + str(count) + "\t")
        x = the_mass_list[i]
        one_comp = the_comp_list[i]
        one_atom = the_atom_list[i]
        the_isp = the_spectra[0][i]
        lost_comp = the_spectra[1][i]
        Z_list = the_spectra[2][i]
        the_isp_01 = transform_the_01(the_spectra[3][i])
        com_score = compared_score(exp_isp_01, the_isp_01)
        if com_score < 1:
            # print(com_score)
            pass_count += 1
            continue
        score, match_list, match_theo_list, score_list, score_mass, match_lost_list, match_z_list \
            = match_isotopic(exp_isp, the_isp, ppm, thresh, min_match_num, lost_comp, Z_list, prob)
        if score > 0:
            dict_match_theo = {}
            for i in range(len(match_list)):
                for j in range(len(match_list[i][0])):
                    mass = match_list[i][0][j]
                    exp_int = match_list[i][1][j]
                    the_int = match_theo_list[i][1][j]
                    if mass not in dict_match_exp.keys():
                        dict_match_exp[mass] = exp_int
                    dict_match_theo[mass] = the_int
            dict_theo_list.append(dict_match_theo)
            all_comp.append(one_comp)
            all_atom.append(one_atom)
            all_comp_mass.append(x)
            all_comp_score.append(score)
            all_comp_match.append(match_list)
            all_score_list.append(score_list)
            all_score_mass.append(score_mass)
            all_comp_lost.append(match_lost_list)
            all_comp_z.append(match_z_list)
        #     print(count, one_comp, x, score)
        #     print(len(dict_theo_list))
        # if count % 100 == 0:
        #     print(count)

    print('passed', pass_count)
    # print("ID\t" + str(count + pass_count) + "\t")
    print('ma1', time() - ttt)
    return all_comp, all_atom, all_comp_mass, all_comp_score, dict_match_exp, \
           dict_theo_list, all_score_list, all_score_mass, all_comp_lost, all_comp_z


def transform_data(dict_exp_match, dict_theo_list, all_comp_score, candidate_sort_index):
    feature = []
    label = []
    mass_list = []
    tmp_comp = np.array([all_comp_score, dict_theo_list]).T
    sort_comp = sorted(tmp_comp, key=lambda x: x[0], reverse=True)

    for mass in dict_exp_match.keys():
        y = dict_exp_match[mass]
        if y <= 0:
            continue
        else:
            x = []
            for i in range(len(dict_theo_list)):
                ind = candidate_sort_index[i]
                if mass in sort_comp[ind][1].keys():
                    tmp_int = sort_comp[ind][1][mass]
                else:
                    tmp_int = 0
                x.append(tmp_int)
            feature.append(x)
            label.append(y)
            mass_list.append(mass)
    feature = np.array(feature)
    label = np.array(label)
    return feature, label, mass_list


def lasso_reg(feature, label, alpha=0.1):
    lasso = Lasso(alpha=alpha, positive=True, max_iter=5000)
    lasso.fit(feature, label)
    y_pred = lasso.predict(feature)
    r2_score_lasso = r2_score(label, y_pred)
    w = lasso.coef_
    print('R2 score:', r2_score_lasso)
    if alpha == 0.5:
        pd.DataFrame(w).to_csv('result/w_0.5.csv', header=None, index=None)
    pl(w)
    return w, r2_score_lasso


def get_label(match_list):
    all_comp = match_list[0]
    all_comp_mass = match_list[2]
    # all_comp_score = match_list[3]
    all_score_list = match_list[6]
    all_score_mass = match_list[7]
    all_comp_lost = match_list[8]
    all_comp_z = match_list[9]

    ## key:理论中性质量，value: 理论结构分子组成
    dict_mass_comp = {}
    ## key:理论中性质量，value: flag 数组， flag[0]: 实验谱匹配m/z,OR 理论中性mass-1, flag[1]: 电荷数，OR 0
    dict_mass_flag = {}
    ## key:理论中性质量，value: 衍生 label 数组 index
    dict_mass_family = {}

    label = []
    for i in range(len(all_comp)):
        tmp_mass = all_comp_mass[i]
        flag = [tmp_mass - 1, 0]
        dict_mass_comp[tmp_mass] = all_comp[i]
        index_list = []
        for j in range(len(all_score_list[i])):
            tmp_mz = all_score_mass[i][j][0]
            tmp_z = all_comp_z[i][j]
            tmp_comp = all_comp[i]
            tmp_lost = all_comp_lost[i][j]
            if np.sum(tmp_lost) == 0:
                flag = [tmp_mz, tmp_z]
            tmp_score = all_score_list[i][j]
            label.append([tmp_mz, tmp_z, tmp_comp, tmp_lost, tmp_score])
            index_list.append(len(label) - 1)
        dict_mass_family[tmp_mass] = index_list
        dict_mass_flag[tmp_mass] = flag
    return label, dict_mass_comp, dict_mass_flag, dict_mass_family


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


def get_n_label(match_list, total_int, candidate_sort_index, n):
    comp = np.array(match_list[0:4] + match_list[5:10]).T
    sorted_comp = sorted(comp, key=lambda x: x[3], reverse=True)
    dict_exp_match = match_list[4]
    ## key:理论中性质量，value: 理论结构分子组成
    dict_mass_comp = {}
    ## key:理论中性质量，value: flag 数组， flag[0]: 实验谱匹配m/z,OR 理论中性mass-1, flag[1]: 电荷数，OR 0
    dict_mass_flag = {}
    ## key:理论中性质量，value: 衍生 label 数组 index
    dict_mass_family = {}
    key_with_order = []
    label = []
    match_int = 0
    for i in dict_exp_match.keys():
        match_int += dict_exp_match[i]
    if len(sorted_comp) <= 0:
        print('No matched candidate!')
        return label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order
    if n < 0 or n > len(candidate_sort_index):
        print('invalid candidate number!')
        return label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order
    score_coff = match_int/total_int
    dict_has_labeled = {}
    max_score = sorted_comp[0][3]
    for i in range(n):
        id = candidate_sort_index[i]
        key_id = id * 10000
        tmp_mass = sorted_comp[id][2]
        tmp_score = np.round(sorted_comp[id][3] * score_coff/max_score, 4)
        flag = [tmp_mass - 1, 0, []]
        new_key = np.round(key_id + tmp_mass, 5)
        dict_mass_comp[new_key] = sorted_comp[id][0]
        key_with_order.append(new_key)
        index_list = []
        all_score_mass = sorted_comp[id][6]
        all_score_list = sorted_comp[id][5]
        max_mass, max_int = get_high_peak(all_score_mass, dict_exp_match)
        is_labeled = False
        max_index = np.array(range(len(max_mass)), dtype=int)
        sort_mass = sorted(np.array([max_mass, max_int, max_index]).T, key=lambda x: x[1], reverse=True)
        k = 0
        while k < len(sort_mass):
            mz = sort_mass[k][0]
            idx = np.int(sort_mass[k][2])
            if mz not in dict_has_labeled.keys():
                z = sorted_comp[id][8][idx]
                lost = sorted_comp[id][7][idx]
                flag = [mz, z, lost, tmp_score]
                is_labeled = True
                dict_has_labeled[mz] = 1
                break
            else:
                k += 1
        if not is_labeled:
            mz = sort_mass[0][0]
            idx = np.int(sort_mass[0][2])
            z = sorted_comp[id][8][idx]
            lost = sorted_comp[id][7][idx]
            flag = [mz, z, lost, tmp_score]
        for j in range(len(max_mass)):
            # tmp_mz = max_mass[j]
            tmp_mz = all_score_mass[j][0]
            order_peak = 1
            for i in (range(len(all_score_mass[j]))):
                mass = all_score_mass[j][i]
                if dict_exp_match[mass] > 0:
                    tmp_mz = mass
                    order_peak = i + 1
                    break
            tmp_score_p = np.round(all_score_list[j], 4)
            tmp_z = sorted_comp[id][8][j]
            tmp_comp = sorted_comp[id][0]
            tmp_lost = sorted_comp[id][7][j]
            label.append([tmp_mz, tmp_z, tmp_comp, tmp_lost, tmp_score, order_peak])
            index_list.append(len(label) - 1)
        dict_mass_family[new_key] = index_list
        dict_mass_flag[new_key] = flag
    return label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order


#  根据额外匹配峰的数目挑选最有信息的结构
def get_most_power_candidate(sorted_comp, has_considered_peak, has_added_candidate, dict_exp):
    max_index = -1
    max_score = -1
    max_peak_num = -1
    for i in range(len(sorted_comp)):
        if i in has_added_candidate:
            continue
        else:
            tmp_score = sorted_comp[i][3]
            tmp_peak_num = 0
            score_mass = sorted_comp[i][-1]
            for j in range(len(score_mass)):
                one_isp_score = score_mass[j]
                for mass in one_isp_score:
                    if (mass not in has_considered_peak) and dict_exp[mass] > 0:
                        tmp_peak_num += 1
            if tmp_peak_num == 0:
                continue
            if tmp_peak_num > max_peak_num:
                max_peak_num = tmp_peak_num
                max_index = i
                max_score = tmp_score
            elif tmp_peak_num == max_peak_num:
                if max_score < tmp_score:
                    max_index = i
                    max_score = tmp_score
    return max_index


# 跟据额外匹配峰的丰度和挑选最有信息的结构
def get_most_power_candidate2(sorted_comp, has_considered_peak, has_added_candidate, dict_exp):
    max_index = -1
    max_score = -1
    max_int_sum = -1
    for i in range(len(sorted_comp)):
        if i in has_added_candidate:
            continue
        else:
            tmp_score = sorted_comp[i][3]
            tmp_int_sum = 0
            score_mass = sorted_comp[i][-1]
            for j in range(len(score_mass)):
                one_isp_score = score_mass[j]
                for mass in one_isp_score:
                    if (mass not in has_considered_peak) and dict_exp[mass] > 0:
                        tmp_int_sum += dict_exp[mass]
            if tmp_int_sum == 0:
                continue
            if tmp_int_sum > max_int_sum:
                max_int_sum = tmp_int_sum
                max_index = i
                max_score = tmp_score
            elif tmp_int_sum == max_int_sum:
                if max_score < tmp_score:
                    max_index = i
                    max_score = tmp_score
    return max_index


# 跟据JS散度值挑选最有信息的结构
def get_most_power_candidate3(sorted_comp, has_considered_peak, has_added_candidate, dict_score, dict_exp):
    max_index = -1
    max_score = -1
    for i in range(len(sorted_comp)):
        if i in has_added_candidate:
            continue
        else:
            dict_tmp = dict_score.copy()
            tmp_score_sum = 0
            score_mass = sorted_comp[i][-1]
            score_list = sorted_comp[i][-2]
            for j in range(len(score_mass)):
                one_isp_mass = score_mass[j]
                first_mass = one_isp_mass[0]
                for mass in one_isp_mass:
                    if dict_exp[mass] > 0:
                        first_mass = mass
                        break
                if first_mass in dict_tmp.keys():
                    tmp_score = dict_tmp[first_mass]
                    if tmp_score < score_list[j]:
                        dict_tmp[first_mass] = score_list[j]
                        tmp_score_sum += 0.0 * (score_list[j] - tmp_score)
                else:
                    dict_tmp[first_mass] = score_list[j]
                    tmp_score_sum += score_list[j]
            if tmp_score_sum == 0:
                continue
            if tmp_score_sum > max_score:
                max_index = i
                max_score = tmp_score_sum

    return max_index


def get_most_power_candidate_score(match_list):
    result_score = []
    s1 = 0
    s2 = 0
    # constant = (math.e*math.e)/(0.1+math.e*math.e)
    constant = 0.99
    r1 = []
    r2 = []
    comp = np.array(match_list[0:4] + match_list[5:8]).T
    sorted_comp = sorted(comp, key=lambda x: x[3], reverse=True)
    n = len(sorted_comp)
    dict_exp = match_list[4]
    has_added_candidate = []
    has_considered_peak = []
    dict_score = {}
    while len(has_added_candidate) < n:
        if len(has_added_candidate) == 0:
            best_index = 0
        else:
            # best_index = get_most_power_candidate(sorted_comp,has_considered_peak,has_considered_peak,dict_exp)
            # best_index = get_most_power_candidate2(sorted_comp,has_considered_peak,has_considered_peak,dict_exp)
            best_index = get_most_power_candidate3(sorted_comp, has_considered_peak, has_considered_peak, dict_score,
                                                   dict_exp)
        if best_index >= 0:
            score_mass = sorted_comp[best_index][-1]
            score_list = sorted_comp[best_index][-2]
            for i in range(len(score_mass)):
                one_isp_mass = score_mass[i]
                for mass in one_isp_mass:
                    if (mass not in has_considered_peak) and dict_exp[mass] > 0:
                        has_considered_peak.append(mass)
                        s1 += dict_exp[mass]
                first_mass = one_isp_mass[0]
                for mass in one_isp_mass:
                    if dict_exp[mass] > 0:
                        first_mass = mass
                        break
                if first_mass in dict_score.keys():
                    tmp_score = dict_score[first_mass]
                    if tmp_score < score_list[i]:
                        dict_score[first_mass] = score_list[i]
                        s2 += 0 * (score_list[i] - tmp_score)
                else:
                    dict_score[first_mass] = score_list[i]
                    s2 += score_list[i]

            has_added_candidate.append(best_index)
            count_c = len(has_added_candidate)
            weight = np.power(constant, count_c)
            r1.append(weight * s1)
            r2.append(weight * s2)
        else:
            # print('Has selected all powerful candidates! ')
            # print('selected/all',len(has_added_candidate),'/',n)
            # break
            count_c = len(has_added_candidate)
            for i in range(n - len(has_added_candidate)):
                weight = np.power(constant, count_c + i)
                r1.append(weight * s1)
                r2.append(weight * s2)
            for i in range(n):
                if i not in has_added_candidate:
                    has_added_candidate.append(i)
            break
    result_score.append(r1)
    result_score.append(r2)
    return result_score, has_added_candidate


def get_accumlate_score(match_list):
    '''

    :param match_list:   匹配返回的结果
    :return: result_score: 累计得分结果
    r1,r2,r3: 三种打分方式的累计得分列表
    s1,s2,s3: 三种打分方式的当前累计得分
    '''
    result_score = []
    r1 = []
    r2 = []
    r3 = []
    s1 = 0
    s2 = 0
    s3 = 0
    has_cosider = []
    dict_exp = match_list[4]
    dict_score = {}
    comp = np.array(match_list[0:4] + match_list[5:8]).T
    sorted_comp = sorted(comp, key=lambda x: x[3], reverse=True)
    for i in range(len(sorted_comp)):
        s1 += sorted_comp[i][3]
        r1.append(s1)
        one_theo_match = sorted_comp[i][-1]
        for x in one_theo_match:
            for mass in x:
                if mass in has_cosider:
                    continue
                else:
                    has_cosider.append(mass)
                    if dict_exp[mass] > 0:
                        s2 += 1
        r2.append(s2)
        score_list = sorted_comp[i][-2]
        score_mass = sorted_comp[i][-1]
        for j in range(len(score_mass)):
            mass = np.round(score_mass[j][0], 4)
            if mass in dict_score.keys():
                tmp_score = dict_score[mass]
                if tmp_score <= score_list[j]:
                    s3 -= tmp_score
                    s3 += score_list[j]
                    dict_score[mass] = score_list[j]
            else:
                s3 += score_list[j]
                dict_score[mass] = score_list[j]
        r3.append(s3)
    result_score.append(r1)
    result_score.append(r2)
    result_score.append(r3)
    return result_score

def get_total_intensity(filter_mz,max_int):
    total_int = 0
    for peak in filter_mz:
        total_int+=peak[1]
    total_int = np.round(total_int/max_int*100,4)
    return total_int

'''
example:
result1=
[[356.123, 321],
 [357.123, 241],
 [358.123, 567],
 [671.234, 6759],
 [672.234, 673],
 [673.234, 534]]

result2=
[[358.123, 567, '-1','3'],
 [671.234, 6759, '-1', '1']
 [495.234, 546, '-1,-2', '1,1']]

result3=
[[356.123, 567, 1, [13, 10, 2, 10, 2], 356.132],
 [671.234, 6759, 1, [16, 20, 3, 17, 4], 671.232]
 ]
 
result4 ={
321.2: [[[321.2, 331.2, 341,2], [10, 11, 12], 10], 
        [[221.2, 321.2],[123, 10], 100]]
}

the_HP: 5个元素的list
the_HP[0] = dict_mass_comp
the_HP[1] = dict_mass_atom
the_HP[2] = all_comp_list
the_HP[3] = all_atom_list
the_HP[4] = all_mass_list
dict_mass_comp ={
351: [[1,4,3,1,2,1,0],[1,3,3,1,2,1,0]]
}
dict_mass_atom ={
351: [[14,12,2,13,2],[14,12,2,13,2]]
}

label：
[[mass, Z, component, lose , score],[mass, Z component, lose , score]]
mass: float, 实验谱匹配峰m/z
Z: int, 电荷数 （正整数）
component: list  肝素组成，7元组
lose: list 各基团丢失数， 4元组
score: float, 结构得分
'''

'''
all_comp : 所有匹配上的结构的基团组成
[[1,2,3,2,2,0,1],[1,3,3,2,1,1,1]]
all_atom : 所有匹配上的结构的元素组成
[[12,11,2,4,2],[14,12,2,10,3]]
all_comp_mass: 所有匹配上的结构的中性质量
[452.1123, 658.4563]
all_comp_score: 所有匹配上的结构得分
[0.3425, 1.4652]
dict_exp_match: 所有结构在实验谱上的匹配上的峰的并集
key: exp谱上的m/z 值
value: m/z 值 对应的exp谱的强度值 (相对值0-100)
{456.3421: 45.3241
 891.2342: 1.2342}
dict_theo_list: 所有结构在理论谱上匹配上的峰的信息
list型, 长度为匹配上的结构数目, 每个元素为一个字典 dict_theo_match, 结构与dict_exp_match 相同
代表每个匹配的候选结构与实验谱匹配上的理论谱的峰的信息
all_score_list: 每个结构中每个匹配同位素峰簇的得分
all_score_mass: 每个结构中 匹配上的实验谱上同位素峰簇mass
'''

"""
def data_process(filter_MZ, max_int, fit_list, the_HP, cand_num=0, ppm=20, bound_th=0.001,
                 bound_intensity=300,
                 delta=None):
    from preprocess import get_filter_MZ_pk, input_candidate, pre_process_data
    # 因为2个文件相互引用了，所以在这里import
    from preprocess import get_filter_MZ_pk, input_candidate

    s = time()
    top_n = 5
    # ppm = 10
    prob = 0.95
    isp_thresh = 0.5
    # the_HP = get_comp_pk('data/HP_20.pk')  # 枚举的理论肝素结构
    # filter_MZ_path = "data/filter_MZ.mgf"
    # spectra_path = "data/938.mgf"
    # spectra_path = "data/0.01_938.mgf"
    the_SP_path = 'data/the_sp_0.1.pk'
    result_path = 'result/'
    # fit_list = get_fit_pk('data/Isotope_dist_fit.pk')
    t0 = time()
    print('read:', t0 - s)
    filter_MZ, max_int = pre_process_data(filter_MZ, max_int, fit_list, ppm, bound_th,
                                          bound_intensity)
    # save_file(filter_MZ, 'data/filter_MZ')
    t1 = time()
    print('preprocess', t1 - t0)
    dict_change_start_list, isotopic_list, dict_list = get_isotopic(filter_MZ, fit_list, ppm, 5, 5, isp_thresh)
    t2 = time()
    print('find isotopic:', t2 - t1)

    the_spectra = get_the_sp(the_SP_path, the_HP[3], 1)
    t7 = time()
    print('generate the sp:', t7 - t2)
    dict_list, delta = align_peak(dict_list, the_HP, the_spectra, top_n, ppm)
    exp_isp = get_exp_isp(dict_list, max_int)
    peaks = filter_MZ
    for i in range(0, len(peaks)):
        peaks[i][0] = np.round(peaks[i][0] + delta, 5)
    filter_MZ = peaks
    all_comp, all_atom, all_comp_mass, all_comp_score, dict_exp_match, \
    dict_theo_list, all_score_list, all_score_mass, all_comp_lost, all_comp_z = \
        calculate_all_hp_score(the_spectra, dict_list, the_HP, max_int, ppm, prob=prob)
    t3 = time()
    print('match time:', t3 - t7)
    match_result = [all_comp, all_atom, all_comp_mass, all_comp_score, dict_exp_match, dict_theo_list, all_score_list,
                    all_score_mass, all_comp_lost, all_comp_z]
    save_pk(match_result, result_path + 'one_match_result0.5.pk')
    result2_score, candidate_sort_index = get_most_power_candidate_score(match_result)
    t4 = time()
    print('pick:', t4 - t3)
    label, dict_mass_comp, dict_mass_flag, dict_mass_family,key_with_order = get_n_label(match_result, candidate_sort_index,
                                                                          len(candidate_sort_index))
    label_save = pd.DataFrame(np.array(label))
    label_save.to_csv('result/label_save.csv', index=None)
    t5 = time()
    print('label:', t5 - t4)

    e = time()
    print('Time cost', e - s)
    max_candi, max_score = 0, 0.0
    for i in range(len(result2_score[1])):
        if result2_score[1][i] > max_score:
            max_candi = i
            max_score = result2_score[1][i]
    label_info = get_n_result(match_result, max_candi + 1)
    return [filter_MZ, max_int, exp_isp, match_result, label_info, len(candidate_sort_index), result2_score[1]]
"""


def get_filter_MZ(origin_MZ, max_int, fit_list, the_HP, ppm=20, bound_th=0.001, bound_intensity=300):
    from preprocess import pre_process_data
    # 因为2个文件相互引用了，所以在这里import

    s = time()
    top_n = 5
    isp_thresh = 0.5
    the_SP_path = 'data/the_sp_0.1.pk'
    t0 = time()
    print('read:', t0 - s)
    filter_MZ, max_int = pre_process_data(origin_MZ, max_int, fit_list, ppm, bound_th,
                                          bound_intensity)
    total_int = get_total_intensity(filter_MZ,max_int)
    # save_file(filter_MZ, 'data/filter_MZ')
    t1 = time()
    print('preprocess', t1 - t0)
    dict_change_start_list, isotopic_list, dict_list = get_isotopic(filter_MZ, fit_list, ppm, 5, 5, isp_thresh)
    t2 = time()
    print('find isotopic:', t2 - t1)

    the_spectra = get_the_sp(the_SP_path, the_HP[3], 1)
    t7 = time()
    t7 = time()
    print('generate the sp:', t7 - t2)
    dict_list, delta = align_peak(dict_list, the_HP, the_spectra, top_n, ppm)
    exp_isp = get_exp_isp(dict_list, max_int)
    peaks = filter_MZ
    for i in range(0, len(peaks)):
        peaks[i][0] = np.round(peaks[i][0] + delta, 5)
    filter_MZ = peaks
    return filter_MZ, max_int, total_int, exp_isp, the_spectra, dict_list


def data_process(exp_isp, max_int, total_int, the_spectra, dict_list, the_HP, ppm=20):
    prob = 0.95
    result_path = 'result/'
    all_comp, all_atom, all_comp_mass, all_comp_score, dict_exp_match, \
    dict_theo_list, all_score_list, all_score_mass, all_comp_lost, all_comp_z = \
        calculate_all_hp_score(the_spectra, dict_list, the_HP, max_int, ppm, prob=prob)
    t3 = time()
    # print('match time:', t3 - t7)
    match_result = [all_comp, all_atom, all_comp_mass, all_comp_score, dict_exp_match, dict_theo_list, all_score_list,
                    all_score_mass, all_comp_lost, all_comp_z]
    save_pk(match_result, result_path + 'one_match_result0.5.pk')
    result2_score, candidate_sort_index = get_most_power_candidate_score(match_result)
    t4 = time()
    print('pick:', t4 - t3)
    label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order = get_n_label(match_result,
                                                                                          total_int,
                                                                                          candidate_sort_index,
                                                                                          len(candidate_sort_index))
    label_save = pd.DataFrame(np.array(label))
    label_save.to_csv('result/label_save.csv', index=None)
    t5 = time()
    print('label:', t5 - t4)

    e = time()
    # print('Time cost', e - s)
    # label_info = get_n_result(match_result, max_candi + 1)
    label_info = [label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order]
    if len(result2_score[1]) > 50:
        label_info = [label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order[:50]]
        return [match_result, label_info, 50, result2_score[1][:50]]
    else:
        return [match_result, label_info, len(candidate_sort_index), result2_score[1]]


def get_n_result(match_result, cand_num):
    result2_score, candidate_sort_index = get_most_power_candidate_score(match_result)
    label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order = get_n_label(match_result,
                                                                                          candidate_sort_index,
                                                                                          cand_num)
    label_info = [label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order]
    return label_info


if __name__ == '__main__':
    # 因为2个文件相互引用了，所以在这里import
    from preprocess import get_filter_MZ_pk, input_candidate, pre_process_data

    s = time()
    ppm = 20
    top_n = 5
    prob = 0.95
    isp_thresh = 0.8
    bound_th = 0.001
    bound_intensity = 300
    the_HP = get_comp_pk('data/HP_20.pk')  # 枚举的理论肝素结构
    filter_MZ_path = "data/filter_MZ.mgf"
    spectra_path = "data/938.mgf"
    # candidate_path = "data/HP_20_mass.csv"
    the_SP_path = 'data/the_sp_0.1.pk'
    result_path = 'result/'
    fit_list = get_fit_pk('data/Isotope_dist_fit.pk')
    # candidates = input_candidate(candidate_path)
    t0 = time()
    print('read:', t0 - s)
    filter_MZ, max_int = get_filter_MZ_pk(filter_MZ_path, spectra_path, ppm)
    # filter_MZ, max_int = get_filter_MZ_pk(filter_MZ_path, spectra_path, the_HP[-1], ppm)
    t1 = time()
    print('preprocess', t1 - t0)
    dict_change_start_list, isotopic_list, dict_list = get_isotopic(filter_MZ, fit_list, ppm, 5, 5,
                                                                    isp_thresh)  # 寻找同位素峰，并返回记录不同价态的字典列表
    t2 = time()
    print('find isotopic:', t2 - t1)
    the_spectra = get_the_sp(the_SP_path, the_HP[3], 1)
    t7 = time()
    print('generate the sp:', t7 - t2)
    dict_list, delta = align_peak(dict_list, the_HP, the_spectra, top_n, ppm)
    all_comp, all_atom, all_comp_mass, all_comp_score, dict_exp_match, \
    dict_theo_list, all_score_list, all_score_mass, all_comp_lost, all_comp_z = \
        calculate_all_hp_score(the_spectra, dict_list, the_HP, max_int, ppm, prob=prob)
    t3 = time()
    print('match time:', t3 - t7)
    match_result = [all_comp, all_atom, all_comp_mass, all_comp_score, dict_exp_match, dict_theo_list, all_score_list,
                    all_score_mass, all_comp_lost, all_comp_z]
    save_pk(match_result, result_path + 'one_match_result0.5.pk')
    # with open(result_path+'one_match_result0.5.pk', 'rb') as f:
    #     match_result= pk.load(f)
    # result_score=get_accumlate_score(match_result)
    result2_score, candidate_sort_index = get_most_power_candidate_score(match_result)
    t4 = time()
    print('pick:', t4 - t3)
    label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order = get_n_label(match_result,
                                                                                          candidate_sort_index,
                                                                                          len(candidate_sort_index))
    label_save = pd.DataFrame(np.array(label))
    label_save.to_csv('result/label_save.csv', index=None)
    t5 = time()
    print('label:', t5 - t4)
    e = time()
    print('Time cost', e - s)

    # result = get_one_mass_comp(520.98219, 2, 20, the_HP)
    # mass = get_component_mass([0,2,2,0,5,1,0])
    # print(mass)
