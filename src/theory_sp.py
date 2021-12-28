import numpy as np
import itertools
import os
import pickle as pk
from collections import Counter as ct
from src.atom_mass import dict_atom
from brainpy import isotopic_variants
from src.utils import save_pk
from src.utils import sort_the


def generate_new_composition(filter_list, max_lost_count=20):
    dict_com = {}
    lost_record = []  # SO3, NH, NHSO3, COO, 丢失数目
    lost_list = [[0, 0, 0, 0, 0],  # 不丢失
                 [0, 0, 0, 3, 1],  # 丢失SO3
                 [1, 0, 1, 0, 0],  # 丢失NH
                 # [2, 2, 0, 1, 0],  # 丢失 COCH2
                 # [1, 0, 1, 3, 1],  # 丢失NHSO3
                 [2, 0, 0, 1, 0],  # 丢失H2O
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


## is_HNa: 是否以脱H加Na带负电, 或者 脱H 加 NH3
## 脱H上限为硫酸根数目加羧基数目
def generate_theory_SP(atom, component, prob=0.5, max_ion=5, is_HNa='H'):
    # all_com, lost_record = generate_new_composition([atom], 2)
    all_atom, all_comp, lost_record = [atom],[component],[[0, 0, 0, 0]]
    # pl(all_com)
    lost_comp = []
    H_Na_count = []
    Z_list = []
    all_mz = []
    # charges = [-1]
    charges = -np.array(range(1,max_ion+1))
    # charges = [-1, -2, -3, -4, -5, -6, -7 ]
    for i in range(len(all_atom)):
        x = all_atom[i]
        y = all_comp[i]
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
            if is_HNa !='H':
                add_mass = 0
                if is_HNa =='HNa':
                    add_mass =dict_atom['Na']-dict_atom['H1']
                elif is_HNa == 'HNH4':
                    add_mass =dict_atom['N14']+4*dict_atom['H1']-dict_atom['H1']
                max_H = y[0]+y[1]+y[4] # 糖醛酸数目加硫酸根数目
                for n in range(-j, max_H+1):
                    diff_mass = 0
                    n_na = n+j # Na的数目是H的数目减电荷数
                    ## for Na
                    # diff_mass = n_na *(dict_atom['Na']-1) # 计算与只脱H的mass的差距
                    ## for NH3
                    # diff_mass = n_na *(dict_atom['N14']+4*dict_atom['H1']-1) # 计算与只脱H的mass的差距
                    diff_mass = n_na * add_mass
                    diff_mass = diff_mass/np.abs(j)
                    tmp_mz = []
                    tmp_int = []
                    for peak in tmp_distribution:
                        tmp_mz.append(np.round(peak.mz+diff_mass, 4))
                        tmp_int.append(peak.intensity * 100 * np.power(prob, sum(losted)))
                        # all_mz.append([np.round(peak.mz,4), peak.intensity])
                    H_Na_count.append([n,n_na])
                    all_mz.append([tmp_mz, np.round(tmp_int, 4)])
                    lost_comp.append(losted)
                    Z_list.append(-j)
            else:
                for peak in tmp_distribution:
                    tmp_mz.append(np.round(peak.mz, 4))
                    tmp_int.append(peak.intensity * 100 * np.power(prob, sum(losted)))
                    # all_mz.append([np.round(peak.mz,4), peak.intensity])
                all_mz.append([tmp_mz, np.round(tmp_int, 4)])
                lost_comp.append(losted)
                H_Na_count.append([-j, 0])
                Z_list.append(-j)
    # sorted_mz = sorted(all_mz, key=lambda mz: mz[0])
    return all_mz, lost_comp, Z_list, H_Na_count


def cal_all_the_sp(all_atom_list, all_comp_list, prob, max_ion=5, is_HNa='H'):
    all_the = []
    all_lost = []
    all_z = []
    all_the_01 = []
    all_HNa = []
    for i in range(len(all_atom_list)):
        one_atom = all_atom_list[i]
        one_comp = all_comp_list[i]
        the_isp, lost_comp, Z_list, H_Na_list = generate_theory_SP(one_atom,one_comp, prob, max_ion, is_HNa=is_HNa)
        the_isp, lost_comp, Z_list, H_Na_list = sort_the(the_isp, lost_comp, Z_list, H_Na_list)
        min_id = 0
        for i in range(len(the_isp)):
            x = the_isp[i][0][0]
            if x >= 100:
                min_id = i
                break
        the_isp = the_isp[min_id:]
        lost_comp = lost_comp[min_id:]
        Z_list = Z_list[min_id:]
        H_Na_list = H_Na_list[min_id:]
        all_the.append(np.float32(the_isp))
        # all_the.append(the_isp)
        all_lost.append(lost_comp)
        all_z.append(Z_list)
        one_01 = change_the_format(the_isp)
        all_the_01.append(one_01)
        all_HNa.append(H_Na_list)
    result = [all_the, all_lost, all_z, all_the_01, all_HNa]
    return result


def get_the_sp(dir, all_atom_list,all_comp_list, prob, max_ion, is_HNa='H'):
    dir = dir + str(max_ion)+is_HNa+'.pk'
    if os.path.exists(dir):
        with open(dir, 'rb') as f:
            data = pk.load(f)
        return data
    else:
        data = cal_all_the_sp(all_atom_list, all_comp_list, prob, max_ion=max_ion, is_HNa=is_HNa)
        save_pk(data, dir)
        return data


def change_the_format(sp, win=0.1):
    index_t = []
    min_th = 100 * np.int(1 / win)
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
    result = np.zeros([1100 * np.int(1 / win)], dtype=bool)
    result[sp] = True
    return result


def change_sp_format(sp, win=0.05):
    result = np.zeros([1100 * np.int(1 / win)], dtype=bool)
    min_th = 100 * np.int(1 / win)
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
