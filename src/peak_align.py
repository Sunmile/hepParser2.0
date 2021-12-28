import numpy as np
import math
import pandas as pd


def get_top_n_isptopic(isptopic_list, n):
    top_n = [] # m/z int z
    top_n_isp = []
    for i in range(len(isptopic_list)):
        z = i + 1
        tmp_dict = isptopic_list[i]
        for mass in tmp_dict.keys():
            tmp_isp = tmp_dict[mass]
            tmp_max_int = np.max(tmp_isp[1])
            tmp_max_id = np.argmax(tmp_isp[1])
            tmp_m_z = tmp_isp[0][tmp_max_id]
            if len(top_n) == 0:
                top_n.append([tmp_m_z, tmp_max_int, z])
                top_n_isp.append(tmp_isp)
            else:
                if [tmp_m_z, tmp_max_int, z] in top_n:
                    continue
                if tmp_max_int < top_n[-1][1]:
                    if len(top_n) < n:
                        top_n.append([tmp_m_z, tmp_max_int, z])
                        top_n_isp.append(tmp_isp)
                    else:
                        continue
                for j in range(1, len(top_n)):
                    if top_n[-j][1] < tmp_max_int < top_n[-(j + 1)][1]:
                        top_n.insert(-j, [tmp_m_z, tmp_max_int, z])
                        top_n_isp.insert(-j, tmp_isp)
                        break
                if tmp_max_int > top_n[0][1]:
                    top_n.insert(0, [tmp_m_z, tmp_max_int, z])
                    top_n_isp.insert(0, tmp_isp)
            if len(top_n) > n:
                top_n.pop()
                top_n_isp.pop()
    max_inten = top_n[0][1]
    return top_n, top_n_isp, max_inten


def change_sp_format1(sp, win=0.1):
    result = np.zeros([1100*np.int(1/win)], dtype=bool)
    min_th = 100*np.int(1/win)
    max_th = 1200*np.int(1/win)
    for one_isp in sp:
        for i in range(len(one_isp[0])):
            mz = one_isp[0][i]
            int = one_isp[1][i]
            if int > 0:
                mass_t = np.int(mz * np.int(1 / win))
                for j in range(11):
                    mass = mass_t + j - 5
                    if min_th <= mass < max_th:
                        result[mass - min_th] = True
    return result


def compared_score1(a, b):
    # return len([1 for i, j in zip(list(a), list(b)) if i==j==True])
    # x = sum((b==True)&(a==True))
    x = np.count_nonzero((b == True) & (a == True))
    return x


def transform_the_011(sp,win=0.1):
    result = np.zeros([1100 * np.int(1 / win)], dtype=bool)
    result[sp]=True
    return result


def sort_exp1(exp_isp):
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

def match_align(exp_isp, the_isp, ppm=10, standard_mz = 1000):
    diff_top_n=np.ones([len(exp_isp)])
    MAX_L = 0
    for item in exp_isp:
        temp = item[0][-1] - item[0][0]
        if temp > MAX_L:
            MAX_L = temp

    for item in the_isp:
        temp = item[0][-1] - item[0][0]
        if temp > MAX_L:
            MAX_L = temp
    i = 0
    j0 = 0
    while (i < len(exp_isp)):
        exp_item = exp_isp[i]
        while (exp_item[0][0] > the_isp[j0][0][0] + MAX_L):
            j0 = j0 + 1
            if j0 == len(the_isp):
                j0 = j0 - 1
                break

        if j0 == len(the_isp) - 1 and exp_item[0][0] > the_isp[j0][0][0] + MAX_L:
            break
        j = j0
        while (j < len(the_isp)):
            the_item = the_isp[j]
            if the_item[0][0] > exp_item[0][0] + MAX_L:
                i += 1
                break
            counter = 0
            tmp_diff = []
            match_str = ''
            win = ppm * the_item[0][0] * 2 * 0.000001
            for aim_mz in the_item[0]:
                match_flag = False
                for k in range(counter, len(exp_item[0])):
                    exp_mz = exp_item[0][k]
                    if abs(standard_mz / aim_mz * (exp_mz - aim_mz)) <= win and exp_item[1][k] != 0:
                        tmp_diff.append(standard_mz / aim_mz * (aim_mz - exp_mz))
                        counter = counter + 1
                        match_str += '1'
                        match_flag = True
                        break
                if not match_flag:
                    match_str += '0'
            if counter < 3 or match_str == '10101':
                j = j + 1
                if j == len(the_isp):
                    i = i + 1
                continue
            diff_top_n[i] = np.mean(tmp_diff)
            i = i + 1
            break
    return diff_top_n


def count_matched(diff, win, delta, top_n, max_int):
    delta = -delta
    count = 0
    sum_diff = 0
    m = max_int
    for i in range(len(diff)):
        min_d = 1
        flag = np.round(np.log10(top_n[i][1] / m * 100 + 1), 4)
        for j in range(len(diff[i])):
            if np.abs(diff[i][j] + delta) < win:
                count = count + flag
                flag = 0
                if np.abs(diff[i][j] + delta) < min_d:
                    min_d = np.abs(diff[i][j] + delta)
        if min_d < 1:
            sum_diff += min_d
    if count == 0:
        sum_diff = 1000
    return count, sum_diff


def get_delta(diff, ppm, top_n, max_int):
    win = ppm * 0.001
    best_delta = 0
    max_count = 0
    min_sum = 1000
    for delta in np.arange(0.1, -0.1, -0.001):
        count, sum_diff = count_matched(diff, win, delta, top_n, max_int)
        if count > max_count:
            max_count = count
            best_delta = delta
            min_sum = sum_diff
        elif count == max_count:
            if sum_diff < min_sum:
                best_delta = delta
                min_sum = sum_diff
    print(max_count)
    return best_delta


def match_top_n(top_n_isp, the_HP, the_spectra, ppm, top_n, max_int):
    top_n_isp = sort_exp1(top_n_isp)
    top_n = sorted(top_n, key=lambda x: x[0])
    top_n_01 = change_sp_format1(top_n_isp)
    all_diff_top_n = []
    pass_count = 0
    count = 0
    tmp_comp_list = the_HP[2]
    for i in range(len(the_spectra[0])):
        one_comp = tmp_comp_list[i]
        # if one_comp == [1,0,1,0,3,0,0]:
        #     a = 1
        the_isp = the_spectra[0][i]
        the_sp_01 = transform_the_011(the_spectra[3][i])
        com_score = compared_score1(top_n_01, the_sp_01)
        if com_score < 1:
            pass_count += 1
            # print('pass',pass_count)
            continue
        diff_top_n = match_align(top_n_isp, the_isp, ppm, 1000)
        # print(count, diff_top_n)
        count += 1
        if np.sum(diff_top_n) < len(diff_top_n):
            all_diff_top_n.append(diff_top_n)
    print('pass', pass_count)
    all_diff_top_n = np.array(all_diff_top_n).T
    # print('diff_matrix',all_diff_top_n)
    diff = []
    for i in range(len(all_diff_top_n)):
        diff.append(np.unique(all_diff_top_n[i]))
    delta = get_delta(diff, ppm, top_n, max_int)
    return np.round(delta, 4)


def get_aligned_isp(isptopic_list, delta):
    result = []
    for i in range(len(isptopic_list)):
        new_dict = {}
        tmp_dict = isptopic_list[i]
        for mass in tmp_dict.keys():
            mz = tmp_dict[mass][0]
            tmp_int = tmp_dict[mass][1]
            new_mz = list(np.array(mz) + delta)
            new_dict[mass + delta] = [np.round(new_mz, 5), tmp_int]
        result.append(new_dict)
    return result


def align_peak(isptopic_list, the_HP, the_spectra, n, ppm):
    top_n, top_n_isp, max_int = get_top_n_isptopic(isptopic_list, n)
    print(top_n)
    delta = match_top_n(top_n_isp, the_HP, the_spectra, ppm, top_n, max_int)
    # delta = 0
    new_isp = get_aligned_isp(isptopic_list,delta)
    print("best delta",delta)
    return new_isp,delta
