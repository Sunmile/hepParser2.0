import numpy as np
from src.utils import get_high_peak


def get_n_label(match_list, total_int, candidate_sort_index, n):
    comp = np.array(match_list[0:4] + match_list[5:11]).T
    sorted_comp = sorted(comp, key=lambda x: x[3], reverse=True)
    dict_exp_match = match_list[4]
    # key:理论中性质量，value: 理论结构分子组成
    dict_mass_comp = {}
    # key:理论中性质量，value: flag 数组， flag[0]: 实验谱匹配m/z,OR 理论中性mass-1, flag[1]: 电荷数，OR 0, flag[2] 脱H加Na
    dict_mass_flag = {}
    # key:理论中性质量，value: 衍生 label 数组 index
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
    score_coff = match_int / total_int
    dict_has_labeled = {}
    max_score = sorted_comp[0][3]
    for i in range(n):
        id = candidate_sort_index[i]
        key_id = id * 10000
        tmp_mass = sorted_comp[id][2]
        # tmp_score = np.round(sorted_comp[id][3] * score_coff / max_score, 4)
        tmp_score = np.round(sorted_comp[id][3],4)
        flag = [tmp_mass - 1, 0,[], [], 0]
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
                HNa = sorted_comp[id][9][idx]
                lost = sorted_comp[id][7][idx]
                flag = [mz, z, HNa, lost, tmp_score]
                is_labeled = True
                dict_has_labeled[mz] = 1
                break
            else:
                k += 1
        if not is_labeled:
            mz = sort_mass[0][0]
            idx = np.int(sort_mass[0][2])
            z = sorted_comp[id][8][idx]
            HNa = sorted_comp[id][9][idx]
            lost = sorted_comp[id][7][idx]
            flag = [mz, z, HNa, lost, tmp_score]
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
            tmp_HNa = sorted_comp[id][9][j]
            tmp_comp = sorted_comp[id][0]
            tmp_lost = sorted_comp[id][7][j]
            label.append([tmp_mz, tmp_z,tmp_HNa, tmp_comp, tmp_lost, tmp_score, order_peak])
            index_list.append(len(label) - 1)
        dict_mass_family[new_key] = index_list
        dict_mass_flag[new_key] = flag
    return label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order


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


# 如果有多个得分相同的结构，会被同时选中
def get_most_power_candidate_score(match_list):
    result_score = []
    s1 = 0
    s2 = 0
    mul_len = 0 # 去除加了好几个得分相同结构时对分数的影响
    # constant = (math.e*math.e)/(0.1+math.e*math.e)
    constant = 0.99
    r1 = []
    r2 = []
    has_added_candidate = []
    has_considered_peak = []
    dict_score = {}
    comp = np.array(match_list[0:4] + match_list[5:8]).T
    sorted_comp = sorted(comp, key=lambda x: x[3], reverse=True)
    n = len(sorted_comp)
    dict_id = {}
    sc_group = [0]
    if n<1:
        print('no matched candidates')
        result_score=[[],[]]
        return result_score, has_added_candidate
    pre_sc = sorted_comp[0][3]
    for i in range(1,n):
        tmp_sc = sorted_comp[i][3]
        if tmp_sc<pre_sc:
            for x in sc_group:
                dict_id[x] = sc_group
            sc_group=[i]
        else:
            sc_group.append(i)
        pre_sc = tmp_sc
    for x in sc_group:
        dict_id[x] = sc_group
    dict_exp = match_list[4]
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

            has_added_candidate.extend(dict_id[best_index])
            mul_len += len(dict_id[best_index])-1
            count_c = len(has_added_candidate)-mul_len
            weight = np.power(constant, count_c)
            for i in range(len(dict_id[best_index])):
                r1.append(weight * s1)
                r2.append(weight * s2)
        else:
            # print('Has selected all powerful candidates! ')
            # print('selected/all',len(has_added_candidate),'/',n)
            # break
            count_c = len(has_added_candidate)-mul_len
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
    """
    :param match_list:   匹配返回的结果
    :return: result_score: 累计得分结果
    r1,r2,r3: 三种打分方式的累计得分列表
    s1,s2,s3: 三种打分方式的当前累计得分
    """
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


def get_n_result(match_result, cand_num):
    result2_score, candidate_sort_index = get_most_power_candidate_score(match_result)
    label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order = get_n_label(match_result,
                                                                                          candidate_sort_index,
                                                                                          cand_num)
    label_info = [label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order]
    return label_info
