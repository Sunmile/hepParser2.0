import numpy as np
import pandas as pd
import copy
from src.atom_mass import dict_atom
import pickle as pk
import matplotlib.pyplot as plt
import seaborn as sns

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

# 根据分子式计算分子质量
def get_molecular_mass(atom_list):
    Mass = 0
    atom_mass = [dict_atom['H1'], dict_atom['C12'], dict_atom['N14'],
                 dict_atom['O16'], dict_atom['S32'], dict_atom['Na']]
    for i in range(len(atom_list)):
        Mass += atom_mass[i] * atom_list[i]
    return np.round(Mass, 5)


def get_component_mass(one_comp):
    comp_atoms = [
        [8, 6, 0, 6, 0],  # 不饱和糖醛酸 HexA
        [10, 6, 0, 7, 0],  # 饱和糖醛酸  GlcA
        [13, 6, 1, 5, 0],  # 葡萄糖胺   GlcN
        [4, 2, 0, 2, 0],  # 乙酰基    Ac
        [2, 0, 0, 4, 1],  # 硫酸基    SO3
        [11, 6, 1, 4, 0],  # 内醚糖    Levoglucosan
        [12, 6, 0, 5, 0]  # 甘露糖(dicp jin)    Man
    ]
    comp_atoms = np.array(comp_atoms)
    tmp_atoms = np.array([2, 0, 0, 1, 0])  # 带一个水
    H2O = np.array([2, 0, 0, 1, 0])  # 任两个基团结合，都脱一个水
    for i in range(len(one_comp)):
        tmp_atoms = tmp_atoms + comp_atoms[i] * one_comp[i] - H2O * one_comp[i]
    mass = get_molecular_mass(tmp_atoms)
    tmp_atoms = tmp_atoms.tolist()
    return np.round(mass, 5), tmp_atoms


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
        # print(count, mass, x)
    return dict_mass_comp, dict_mass_atom, all_list, all_atoms_list, all_mass_list


def draw(dir, r2, exp_sp, fig_dir):
    max_int = max(exp_sp[:, 1])
    an_01 = pd.read_excel(dir)
    plt.figure(figsize=(7, 4), dpi=250)
    lw1 = 0.5
    lw2 = 0.5
    for i in range(an_01.shape[0]):
        x = an_01.loc[i, 'mz']
        y1 = an_01.loc[i, 'exp_int']
        y2 = an_01.loc[i, 'the_int']
        w = an_01.loc[i, 'weight']
        w = w[1:-1]
        w_list = w.split(', ')
        if y2>0:
            w_list = [np.float(x) for x in w_list]
        if len(w_list) ==1:
            plt.plot([x, x], [0, -y2], lw=lw2, c='b', zorder=2)
        elif len(w_list) == 2:
            tmp_y1 = y2 * (w_list[0]/sum(w_list))
            tmp_y2 = y2 * (w_list[1]/sum(w_list))
            plt.plot([x, x], [0, -tmp_y1], lw=lw2, c='r', zorder=2)
            plt.plot([x, x], [-tmp_y1, -y2], lw=lw2, c='b', zorder=2)
        elif len(w_list) == 3:
            tmp_y1 = y2 * (w_list[0]/sum(w_list))
            tmp_y2 = y2 * (w_list[1]/sum(w_list))
            tmp_y3 = y2 * (w_list[2]/sum(w_list))
            plt.plot([x, x], [0, -tmp_y1], lw=lw2, c='r', zorder=2)
            plt.plot([x, x], [-tmp_y1, -tmp_y1-tmp_y2], lw=lw2, c='b', zorder=2)
            plt.plot([x, x], [-tmp_y1-tmp_y2, -y2], lw=lw2, c='#ffe66d', zorder=2)
        else:
            print('There are ',str(len(w_list)), ' components can explain peak ',str(x))
            plt.plot([x, x], [0, -y2], lw=lw2, c='b', zorder=2)
    for i in range(len(exp_sp)):
        x = exp_sp[i, 0]
        y = exp_sp[i, 1]/max_int * 100
        plt.plot([x, x], [0, y], lw=lw1, c='grey', zorder=1)
    sns.despine()
    plt.plot([0, 0], [0, 0], lw=lw1, c='grey', label='Exp_sp')
    plt.plot([0, 0], [0, 0], lw=lw2, c='b', label='The_sp')
    plt.plot([200, 1200], [0, 0], lw=0.5, c='k')

    plt.xlim(200, 1200)
    plt.ylim(-100, 100)
    plt.xticks(fontsize=9)
    plt.yticks(fontsize=9)
    plt.xlabel('m/z', fontsize=10)
    plt.ylabel('Relative Intensity',fontsize=10)
    plt.legend(frameon=False, fontsize=8)
    plt.title('R square=' + str(np.round(r2, 5)),fontsize=10)
    plt.savefig(fig_dir)
    plt.show()


if __name__ == '__main__':
    dir = './../data/annotate_BEH5_0.1.csv'
    draw(dir)
    # y = [79.12205786,70.10151105,42.24685991,27.98468934,17.22269964,11.28688745]
    # y2 = [78.9214,70.2802,42.5167,28.3007,18.071,11.6137]
    # x = [0,1,2,3,4,5]

    # y = [88.03543732,99.03938312,20.61357166,40.83050654,8.621136652,17.72581981,3.629045612]
    # y2 = [87.4233,100,21.0062,41.6055,8.778,18.9312,3.5018]
    # x = [0,1,2,3,4,5,6]
    # plt.figure(figsize=(6, 3), dpi=250)
    # for i in range(len(y)):
    #     plt.plot([x[i], x[i]],[0, -y[i]], color='b')
    #     plt.plot([x[i], x[i]],[0, y2[i]], color='grey')
    #
    # plt.plot([-1,7],[0,0],color='k')
    # sns.despine()
    # # plt.plot(x,y)
    # plt.show()
