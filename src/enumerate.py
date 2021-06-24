import numpy as np
import pandas as pd
import copy
import os
from src.utils import *
from src.atom_mass import *


# 考虑所有可能的肝素组成，返回所有肝素组成式
def enumerate_component(n=20):
    all_list = []  # bubaohe, baohe, putaotangan, yixian, liusuan, neimi, ganlu
    count = 0
    num_bu = [0, 1]
    for x in num_bu:
        for i in range(0, n):
            if x == 0:
                if i == 0:
                    num_bao = [i, i + 1]
                else:
                    num_bao = [i - 1, i, i + 1]
            else:
                num_bao = [i, i + 1]
            for j in num_bao:
                if j == i + 1:
                    num_nei = [0]
                    num_gan = [0]
                elif j == i - 1:
                    num_nei = [0, 1]
                    num_gan = [0, 1]
                else:
                    if x == 0:
                        num_nei = [0]
                        num_gan = [0]
                    else:
                        num_nei = [0, 1]
                        num_gan = [0, 1]
                for y in num_nei:
                    for z in num_gan:
                        if y == 1 and z == 1:
                            continue
                        for m in range(j):
                            for p in range(3 * j - m + i + x + y + z):
                                tmp_comp = [x, i, j, m, p, y, z]
                                all_list.append(tmp_comp)
                                count += 1
                                print(count, tmp_comp)

    print('All components counts', len(all_list))
    return all_list


# 根据DP值枚举可能的结构
def enumerate_enoxaparin(dp):
    dict_dp_comp = {} # dp和组分的字典
    dict_l_dp_comp = {} # 单数dp 保留左侧 的字典
    dict_r_dp_comp = {} # 单数dp 保留右侧 的字典

    dp2 = [
        [1, 0, 1, 0, 0, 0, 0],
        [1, 0, 1, 0, 1, 0, 0],
        [1, 0, 1, 0, 2, 0, 0],
        [1, 0, 1, 0, 3, 0, 0],
        [1, 0, 1, 0, 4, 0, 0],
        [1, 0, 1, 1, 0, 0, 0],
        [1, 0, 1, 1, 1, 0, 0],
        [1, 0, 1, 1, 2, 0, 0],
        [1, 0, 1, 1, 3, 0, 0],
        [1, 0, 0, 0, 1, 1, 0],
        [1, 0, 0, 0, 2, 1, 0],
    ]  # dp为2时可能 list
    dict_dp_comp[2]=dp2
    dp1l = [
        [1, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 0, 1, 0, 0],
    ]  # dp为1且保留左侧
    dp1r = [
        [0, 0, 0, 0, 1, 1, 0],
        [0, 0, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 1, 0, 0],
        [0, 0, 1, 0, 2, 0, 0],
        [0, 0, 1, 0, 3, 0, 0],
        [0, 0, 1, 1, 0, 0, 0],
        [0, 0, 1, 1, 1, 0, 0],
        [0, 0, 1, 1, 2, 0, 0]
    ]  # dp为1且保留右侧
    dp_cell =[
        [0, 1, 1, 0, 0, 0, 0],
        [0, 1, 1, 0, 1, 0, 0],
        [0, 1, 1, 0, 2, 0, 0],
        [0, 1, 1, 0, 3, 0, 0],
        [0, 1, 1, 0, 4, 0, 0],
        [0, 1, 1, 1, 0, 0, 0],
        [0, 1, 1, 1, 1, 0, 0],
        [0, 1, 1, 1, 2, 0, 0],
        [0, 1, 1, 1, 3, 0, 0]
    ]  # 2糖 单元可能的结构 list

    dp_cell_1 =[
        [0, 1, 1, 0, 0, 0, 0],
        [0, 1, 1, 0, 1, 0, 0],
        [0, 1, 1, 0, 2, 0, 0],
        [0, 1, 1, 0, 3, 0, 0],
        [0, 1, 1, 1, 0, 0, 0],
        [0, 1, 1, 1, 1, 0, 0],
        [0, 1, 1, 1, 2, 0, 0]
    ]  # 第一个2糖单元需特殊考虑
    dict_l_dp_comp[1]=dp1l
    dict_r_dp_comp[1]=dp1r
    dict_dp_comp[1]=dp1l+dp1r
    i = 2
    while i <= dp:
        i+=1
        tmp= dict_dp_comp[i-2]
        result_list = []
        if i <=4:
            for x in tmp:
                for y in dp_cell_1:
                    new = np.array(x) + np.array(y)
                    new = new.tolist()
                    if new not in result_list:
                        result_list.append(new)
        else:
            for x in tmp:
                for y in dp_cell:
                    new = np.array(x)+np.array(y)
                    new = new.tolist()
                    if new not in result_list:
                        result_list.append(new)
        dict_dp_comp[i] = result_list
    for i in range(dp):
        tmp_res = dict_dp_comp[i+1]
        res_df = pd.DataFrame(tmp_res)
        # res_df.to_csv('./data/theory_structure/enoxaparin_'+str(i+1)+'.csv',header=None, index=None)
    return dict_dp_comp


def enumerate_dalteparin(dp):
    dict_dp_comp = {}  # dp和组分的字典
    dict_l_dp_comp = {}  # 单数dp 保留左侧 的字典
    dict_r_dp_comp = {}  # 单数dp 保留右侧 的字典

    dp2 = [
        [0, 1, 0, 0, 0, 0, 1],
        [0, 1, 0, 0, 1, 0, 1],
        [0, 1, 0, 0, 2, 0, 1],
    ]   # dp为2时可能 list
    dict_dp_comp[2] = dp2
    dp1l = [
        [0, 1, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 1, 0, 0],
    ]   # dp为1且保留左侧
    dp1r = [
        [0, 0, 0, 0, 0, 0, 1],
        [0, 0, 0, 0, 1, 0, 1],
    ]   # dp为1且保留右侧
    dp_cell = [
        [0, 1, 1, 0, 0, 0, 0],
        [0, 1, 1, 0, 1, 0, 0],
        [0, 1, 1, 0, 2, 0, 0],
        [0, 1, 1, 0, 3, 0, 0],
        [0, 1, 1, 0, 4, 0, 0],
        [0, 1, 1, 1, 0, 0, 0],
        [0, 1, 1, 1, 1, 0, 0],
        [0, 1, 1, 1, 2, 0, 0],
        [0, 1, 1, 1, 3, 0, 0]
    ]   # 2糖 单元可能的结构 list
    dict_l_dp_comp[1] = dp1l
    dict_r_dp_comp[1] = dp1r
    dict_dp_comp[1] = dp1l + dp1r
    i = 2
    while i < dp:
        i += 1
        tmp = dict_dp_comp[i - 2]
        result_list = []
        for x in tmp:
            for y in dp_cell:
                new = np.array(x) + np.array(y)
                new = new.tolist()
                if new not in result_list:
                    result_list.append(new)
        dict_dp_comp[i] = result_list
    for i in range(dp):
        tmp_res = dict_dp_comp[i+1]
        res_df = pd.DataFrame(tmp_res)
        # res_df.to_csv('./data/theory_structure/dalteparin_'+str(i+1)+'.csv',header=None, index=None)
    return dict_dp_comp


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
        print(count, mass, x)
    return dict_mass_comp, dict_mass_atom, all_list, all_atoms_list, all_mass_list


def get_comp_pk(dir,dp):
    dir = dir + str(dp)+'.pk'
    if os.path.exists(dir):
        with open(dir, 'rb') as f:
            data = pk.load(f)
        return data
    else:
        # all_list = enumerate_component(20)  # 枚举所有可能肝素分子的构成
        # dict_dp_list = enumerate_dalteparin(dp)
        dict_dp_list = enumerate_enoxaparin(dp)
        all_list = dict_dp_list[dp]
        # all_list = filter_component(all_list, 100, 5000)
        dict_mass_comp, dict_mass_atom, all_comp_list, all_atoms_list, all_mass_list = transform_component_to_atom(
            all_list)
        result = [dict_mass_comp, dict_mass_atom, all_comp_list, all_atoms_list, all_mass_list]
        save_pk(result, dir)
        return result

