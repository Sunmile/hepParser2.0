import numpy as np
import pandas as pd
import copy
import itertools
import pickle as pk
from src.enumerate import *
from pprint import pprint as pl
from sklearn.metrics import r2_score
from sklearn.linear_model import Lasso


def get_distance_list(MZ_list):
    distance_list = []
    for j in range(1, len(MZ_list)):
        distance_list.append(MZ_list[j][0] - MZ_list[j - 1][0])
    return distance_list


def generate_isotope_peak(M, Z=1, isotopic_num=4):
    iso_list = []
    for j in range(isotopic_num):
        iso_list.append(M + j / float(Z))
    return iso_list


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


def get_one_mass_comp(mass, z, ppm, the_HP):
    atom_list = the_HP[1]
    comp_mass = mass * z + z
    win = comp_mass * ppm * 0.000001
    result = []
    for m in atom_list.keys():
        one_atom = atom_list[m][0]
        # all_list, lost_record = generate_new_composition([one_atom], 5)
        all_list, lost_record = [one_atom], [[0, 0, 0, 0]]
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


# def get_label(match_list):
#     all_comp = match_list[0]
#     all_comp_mass = match_list[2]
#     # all_comp_score = match_list[3]
#     all_score_list = match_list[6]
#     all_score_mass = match_list[7]
#     all_comp_lost = match_list[8]
#     all_comp_z = match_list[9]
#
#     ## key:理论中性质量，value: 理论结构分子组成
#     dict_mass_comp = {}
#     ## key:理论中性质量，value: flag 数组， flag[0]: 实验谱匹配m/z,OR 理论中性mass-1, flag[1]: 电荷数，OR 0
#     dict_mass_flag = {}
#     ## key:理论中性质量，value: 衍生 label 数组 index
#     dict_mass_family = {}
#
#     label = []
#     for i in range(len(all_comp)):
#         tmp_mass = all_comp_mass[i]
#         flag = [tmp_mass - 1, 0]
#         dict_mass_comp[tmp_mass] = all_comp[i]
#         index_list = []
#         for j in range(len(all_score_list[i])):
#             tmp_mz = all_score_mass[i][j][0]
#             tmp_z = all_comp_z[i][j]
#             tmp_comp = all_comp[i]
#             tmp_lost = all_comp_lost[i][j]
#             if np.sum(tmp_lost) == 0:
#                 flag = [tmp_mz, tmp_z]
#             tmp_score = all_score_list[i][j]
#             label.append([tmp_mz, tmp_z, tmp_comp, tmp_lost, tmp_score])
#             index_list.append(len(label) - 1)
#         dict_mass_family[tmp_mass] = index_list
#         dict_mass_flag[tmp_mass] = flag
#     return label, dict_mass_comp, dict_mass_flag, dict_mass_family


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

def save_filter_MZ(data, dir):
    with open(dir, 'w', newline='') as f:
        f.write("BEGIN IONS\nTITLE=rawdata\nRTINSECONDS=3.0\nCHARGE=1\nmsLevel=1\n")
        for i in range(len(data)):
            f.writelines(str(data[i][0]) + ' ' + str(data[i][1]))
            f.writelines('\n')
        f.write("END IONS\n")


# data_process_fix中有这个函数，一模一样，都用的and
def filter_spectra(MZ_int, th, max_int, min_int=300):
    filter_MZ = []
    for i in range(len(MZ_int)):
        if MZ_int[i][1] > min_int and MZ_int[i][2] > th:
            filter_MZ.append(MZ_int[i])
    return filter_MZ



def filter_spectra_noise(MZ_int, min_int=300):
    filter_MZ = []
    for i in range(len(MZ_int)):
        if MZ_int[i][1] > min_int:
            filter_MZ.append(MZ_int[i])
    return filter_MZ


# 统计匹配数目
def match_peak_ppm(peaks, max_int, candidates, fit_list, ppm):
    """
    匹配：滑动窗口，以候选结构为基准，在左右tolerace范围内匹配，只要实验峰被匹配上就行，不处理一个实验峰匹配上多个理论结构这种情况

    exp: n*2*5维的实验谱，且5维的mz_list是近似于等间隔的中性质量
    peaks: n * 3 [[m/z,absolute_intensity,relative_intensity]]
    candidates: n [mass1,...]
    fit_list: Isotope_dist_fit.pk
    返回 匹配数目
    """
    dict_score_list, dict_list, dict_the_iso_list = get_isotopic(peaks, fit_list, ppm, 5, 5, 0)
    exp_isp = get_exp_isp(dict_list, max_int)
    exp_peaks = []  # 将exp_isp中每个元素的前3个峰合并起来，方便后续使用滑动窗口法
    num_matched_dict = {}

    for i in range(0, len(exp_isp)):
        mz_list, inten_list = exp_isp[i]
        z = int(1 / (mz_list[1] - mz_list[0]) + 0.5)  # 四舍五入
        for j in range(0, 3):
            if inten_list[j] != 0:
                exp_peaks.append([mz_list[j] * z + z, i])  # tmpmass = mass * z + z

    def cmp(peak1, peak2):
        if peak1[0] < peak2[0]:
            return 1
        elif peak1[0] > peak2[0]:
            return -1
        else:
            return 0

    i = 0
    sorted(exp_peaks, key=cmp_to_key(cmp), reverse=1)
    for j in range(0, len(exp_peaks)):
        mass, id = exp_peaks[j]
        while i < len(candidates):
            tolerance = candidates[i] * ppm * ppm_value
            if mass - candidates[i] > tolerance:
                i += 1
            elif candidates[i] - mass <= tolerance:
                if id not in num_matched_dict.keys():
                    num_matched_dict[id] = 1
                i += 1
            else:
                break
    # for key in num_matched_dict.keys():
    #     if exp_isp[key][1][0] >= 0.2 or exp_isp[key][1][1] >= 0.2 or exp_isp[key][1][2] >= 0.2 \
    #             or exp_isp[key][1][3] >= 0.2 or exp_isp[key][1][4] >= 0.2:
    #         print(exp_isp[key])
    return len(num_matched_dict)


# 网格搜索
def align_grid_search(peaks, max_int, candidates, fit_list, window, steps=10, ppm=10):
    # print("-----------grid search start-----------")
    """
    以匹配数目 为目标函数，对差值进行网格搜索,window是网格搜索的范围
    每次在当前最大值的左右window范围内寻找下一个最大值，window_next = window * 2/ steps

    返回 谱整体的偏移量
    peaks:  n*2*5
    """
    mz_delta_cal, num_delta_max = 0, 0
    new_peaks = [peak.copy() for peak in peaks]
    num_matched_delta, num_matched_max, last_matched_max = 0, 0, -1
    num_iter = 0
    while num_matched_max > last_matched_max and num_iter < 3:
        # while num_matched_max > last_matched_max:
        last_matched_max = num_matched_max
        for i in range(0, steps + 1):
            mz_delta = - window / 2 + window * i / steps
            for j in range(0, len(peaks)):
                new_peaks[j][0] = peaks[j][0] + mz_delta_cal + mz_delta
            num_matched = match_peak_ppm(new_peaks, max_int, candidates, fit_list, ppm)
            # print(mz_delta_cal + mz_delta, num_matched)
            if num_matched > num_matched_max:
                num_delta_max = mz_delta_cal + mz_delta
                num_matched_max = num_matched
        mz_delta_cal = num_delta_max
        window = window * 2 / steps
        num_iter += 1
    # print("-----------grid_search finished. delta=", round(num_delta_max, 6), "num_matched=", num_matched_max)
    return np.round(mz_delta_cal, 5)


# 1.峰合并
def merge_peaks_old(peaks):
    """
    对峰分段 + 求质心,前开后闭(区分分割点算在前一个峰里)，这里采用先求出相对强度，再求质心，也可以 先求质心，再计算相对强度
    :param peaks:
    :return:
    """

    max_int = 0
    new_peaks = []
    mz_abs_cal, intensity_abs_cal, count = peaks[0][0] * peaks[0][1], peaks[0][1], 0
    if (peaks[0][1] > 0.0):
        count = 1
    is_up = True
    for i in range(1, len(peaks)):
        if is_up:  # 上升阶段
            mz_abs_cal += peaks[i][0] * peaks[i][1]
            intensity_abs_cal += peaks[i][1]
            if peaks[i][1] != 0.0:
                count += 1
            if peaks[i][1] < peaks[i - 1][1]:  # 上升的拐点
                is_up = False
        else:
            if peaks[i][1] <= peaks[i - 1][1]:  # 下降阶段
                mz_abs_cal += peaks[i][0] * peaks[i][1]
                intensity_abs_cal += peaks[i][1]
                if (peaks[i][1] != 0.0):
                    count += 1
            else:  # 下降的拐点

                # m/z求质心,绝对intensity累加，相对intensity求均值
                # 若此处谷底低于峰的20%以下，截断
                # 若此处高于峰的20%以上，认为还在同一个峰里
                new_peaks.append(
                    [np.round(mz_abs_cal / intensity_abs_cal, 5), intensity_abs_cal])

                max_int = max(max_int, intensity_abs_cal)
                mz_abs_cal, intensity_abs_cal = peaks[i][0] * peaks[i][1], peaks[i][1]
                count = 1
                is_up = True
    return new_peaks, max_int


def merge_peaks(peaks, max_int=0, threshold=0.2):
    """
    对峰分段 + 求质心,前开后闭(区分分割点算在前一个峰里)，这里采用先求出相对强度，再求质心，也可以 先求质心，再计算相对强度

    峰的截断：降低：峰底高度 <= 左侧峰高 *  α
            升高：峰底高度 <= 右侧峰高 *  α
            同时满足降低和升高高度差阈值时，截断峰
    相对只进行波形合并的merge_peaks_old()，峰的数目减少 30%
    """

    threshold = 0.2
    max_int = 0
    new_peaks = []
    count, last_count = 0, 0
    last_peak_top, last_peak_bottom = 0, 0
    mz_abs_cal, last_mz_abs_cal = peaks[0][0] * peaks[0][1], 0
    intensity_abs_cal, last_intensity_abs_cal = peaks[0][1], 0
    if peaks[0][1] > 0.0:
        count = 1
    is_up = True
    for i in range(1, len(peaks)):
        if peaks[i][0]>=348.5:
            aa = 0
        if is_up:  # 上升阶段
            mz_abs_cal += peaks[i][0] * peaks[i][1]
            intensity_abs_cal += peaks[i][1]
            if peaks[i][1] != 0.0:
                count += 1
            if peaks[i][1] < peaks[i - 1][1]:  # 上升的拐点
                # 判断左侧的波谷是否要截断成峰
                if last_peak_bottom < last_peak_top * threshold and last_peak_bottom < peaks[i - 1][1] * threshold:
                    new_peaks.append([np.round(last_mz_abs_cal / last_intensity_abs_cal, 5),
                                      last_intensity_abs_cal])
                    max_int = max(max_int, last_intensity_abs_cal)
                else:
                    count += last_count
                    mz_abs_cal += last_mz_abs_cal
                    intensity_abs_cal += last_intensity_abs_cal
                is_up = False
                last_peak_top = peaks[i - 1][1]
        else:
            if peaks[i][1] <= peaks[i - 1][1]:  # 下降阶段
                mz_abs_cal += peaks[i][0] * peaks[i][1]
                intensity_abs_cal += peaks[i][1]
                if peaks[i][1] != 0.0:
                    count += 1
            else:  # 下降的拐点
                last_peak_bottom = peaks[i - 1][1]
                last_count = count
                last_mz_abs_cal = mz_abs_cal
                last_intensity_abs_cal = intensity_abs_cal

                is_up = True
                count = 1
                mz_abs_cal = peaks[i][0] * peaks[i][1]
                intensity_abs_cal = peaks[i][1]
    # 少最后一个峰
    # new_peaks = peaks
    for i in range(len(new_peaks)):
        new_peaks[i].append(new_peaks[i][1] / max_int)
    return new_peaks, max_int



def merge_peaks_with_original_mz(peaks, max_int=0, threshold=0.2):
    """
    对峰分段 + 求质心,前开后闭(区分分割点算在前一个峰里)，这里采用先求出相对强度，再求质心，也可以 先求质心，再计算相对强度

    峰的截断：降低：峰底高度 <= 左侧峰高 *  α
            升高：峰底高度 <= 右侧峰高 *  α
            同时满足降低和升高高度差阈值时，截断峰
    相对只进行波形合并的merge_peaks_old()，峰的数目减少 30%
    """
    #

    new_peaks = []
    count, last_count = 0, 0
    high_peak_int, high_peak_mz = 0, 0
    last_peak_top, last_peak_bottom = 0, 0
    mz_abs_cal, last_mz_abs_cal = peaks[0][0] * peaks[0][1], 0
    intensity_abs_cal, last_intensity_abs_cal = peaks[0][1], 0
    if peaks[0][1] > 0.0:
        count = 1
    is_up = True
    for i in range(1, len(peaks)):
        if high_peak_int < peaks[i][1]:
            high_peak_mz = peaks[i][0]
            high_peak_int = peaks[i][1]
        if is_up:  # 上升阶段
            mz_abs_cal += peaks[i][0] * peaks[i][1]
            intensity_abs_cal += peaks[i][1]
            if peaks[i][1] != 0.0:
                count += 1
            if peaks[i][1] < peaks[i - 1][1]:  # 上升的拐点
                # 判断左侧的波谷是否要截断成峰
                if last_peak_bottom < last_peak_top * threshold and last_peak_bottom < peaks[i - 1][1] * threshold:
                    new_peaks.append([high_peak_mz, last_intensity_abs_cal])
                    max_int = max(max_int, last_intensity_abs_cal)
                    high_peak_int = 0
                else:
                    count += last_count
                    mz_abs_cal += last_mz_abs_cal
                    intensity_abs_cal += last_intensity_abs_cal
                is_up = False
                last_peak_top = peaks[i - 1][1]
        else:
            if peaks[i][1] <= peaks[i - 1][1]:  # 下降阶段
                mz_abs_cal += peaks[i][0] * peaks[i][1]
                intensity_abs_cal += peaks[i][1]
                if peaks[i][1] != 0.0:
                    count += 1
            else:  # 下降的拐点
                # m/z求质心,绝对intensity累加，相对intensity求均值
                # new_peaks.append([round(mz_abs_cal / intensity_abs_cal, 5), intensity_abs_cal])
                last_peak_bottom = peaks[i - 1][1]
                last_count = count
                last_mz_abs_cal = mz_abs_cal
                last_intensity_abs_cal = intensity_abs_cal

                is_up = True
                count = 1
                mz_abs_cal = peaks[i][0] * peaks[i][1]
                intensity_abs_cal = peaks[i][1]
    # 少最后一个峰
    # new_peaks = peaks
    for i in range(len(new_peaks)):
        new_peaks[i].append(new_peaks[i][1] / max_int)
    return new_peaks, max_int


# 2.校准谱
def align_peaks(peaks, max_int, candidates, ppm, delta):
    """
    peaks：3维：m/z,absolute_intensity,relative_intensity
    exp是n*2*5维的实验谱，且5维的mz_list是近似于等间隔的中性质量
    匹配的时候用exp_isp 最后移的时候移filter_MZ, 然后重新求exp_isp 再匹配
    返回是和peaks相同的格式
    """
    fit_list = get_fit_pk('../data/Isotope_dist_fit.pk')  # 4.寻找同位素峰，并返回记录不同价态的字典列表
    # if delta is None:
    #     delta = align_grid_search(peaks, max_int, candidates, fit_list, window=0.1, steps=10, ppm=ppm)
    # for i in range(0, len(peaks)):
    #     peaks[i][0] = np.round(peaks[i][0] + delta, 5)
    # delta = align_grid_search(peaks, max_int, candidates, fit_list, window=0.1, steps=10, ppm=ppm)
    print('delta',delta)
    for i in range(0, len(peaks)):
        peaks[i][0] = np.round(peaks[i][0] + delta, 5)
    return peaks, delta
    # +1.15036
    # +0.0535

