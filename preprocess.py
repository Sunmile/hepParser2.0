# from data_process import *
from functools import cmp_to_key
from Hp_opt import *
import numpy as np
import os
import time
import csv

ppm_value = 0.000001


def input_candidate(candidate_path):
    candidates = []
    with open(candidate_path, 'r') as f:
        reader = csv.reader(f)
        for line in reader:
            candidates.append(float(line[0]))

    candidates = list(set(candidates))
    candidates.sort()
    return candidates


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


def format_mz(peaks):
    filter_MZ = []
    for i in range(len(peaks)):
        if peaks[i][0] >= 350:
            filter_MZ.append([np.round(peaks[i][0], 4), peaks[i][1]])
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
    dict_change_start_list, isotopic_list, dict_list = get_isotopic(peaks, fit_list, ppm, 5, 5, 0)
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


def merge_peaks(peaks):
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
    fit_list = get_fit_pk('data/Isotope_dist_fit.pk')  # 4.寻找同位素峰，并返回记录不同价态的字典列表
    # if delta is None:
    #     delta = align_grid_search(peaks, max_int, candidates, fit_list, window=0.1, steps=10, ppm=ppm)
    # for i in range(0, len(peaks)):
    #     peaks[i][0] = np.round(peaks[i][0] + delta, 5)
    # delta = align_grid_search(peaks, max_int, candidates, fit_list, window=0.1, steps=10, ppm=ppm)
    print('delta',delta)
    for i in range(0, len(peaks)):
        peaks[i][0] = np.round(peaks[i][0] + delta+0.0535, 5)
    return peaks, delta+0.0565
           # +1.15036
           # +0.0535


# 预处理的整个流程：峰合并、峰过滤、校准谱
def pre_process_data(filter_MZ, max_int, fit_list, ppm=10, bound_th=0.001, bound_intensity=300):
    '''
    # baseline
    filter_MZ = filter_spectra(filter_MZ, 0.01, max_int, 100)  # 过滤谱文件，去除强度小于100且相对丰度低于0.01的峰
    filter_MZ = local_denoise(filter_MZ, 0.1, 0.2)  # 针对局部进行峰过滤，只保留局部最高峰
    num_matched = match_peak_ppm(filter_MZ, max_int, candidates, fit_list, ppm)
    print("baseline:", "max_int=" + str(max_int), "peak_num=" + str(len(filter_MZ)), "matched_num=" + str(num_matched))
    '''

    filter_MZ = format_mz(filter_MZ)
    # 原始峰
    # num_matched = match_peak_ppm(filter_MZ, max_int, candidates, fit_list, ppm)
    # print("origin:", "peak_num=" + str(len(filter_MZ)), "matched_num=" + str(num_matched))

    # 峰合并
    filter_MZ, max_int = merge_peaks(filter_MZ)  # 2.针对局部进行峰过滤，只保留局部最高峰
    # num_matched = match_peak_ppm(filter_MZ, max_int, candidates, fit_list, ppm)
    # print("merge:", "max_int=" + str(max_int), "peak_num=" + str(len(filter_MZ)), "matched_num=" + str(num_matched))

    # 过滤
    filter_MZ = filter_spectra(filter_MZ, bound_th, max_int, bound_intensity)  # 3.过滤谱文件，去除强度小于300且相对丰度低于0.001的峰
    # num_matched = match_peak_ppm(filter_MZ, max_int, candidates, fit_list, ppm)
    # print("filter:", "th = 0.001", "min_int = 300 ", "peak_num=" + str(len(filter_MZ)),
    #       "matched_num=" + str(num_matched))

    # 校准
    # filter_MZ, delta = align_peaks(filter_MZ, max_int, candidates, ppm, delta)  # 4.校准谱
    # num_matched = match_peak_ppm(filter_MZ, max_int, candidates, fit_list, ppm)
    # print("align:", "peak_num=" + str(len(filter_MZ)), "matched_num=" + str(num_matched))
    return filter_MZ, max_int


def get_filter_MZ_pk(filter_MZ_path, spectra_path, ppm=10, bound_th=0.001, bound_intensity=300, delta=None):
    if os.path.exists(filter_MZ_path):
        filter_MZ, max_int = read_spectra(filter_MZ_path)
    else:
        filter_MZ, max_int = read_spectra(spectra_path)
        fit_list = get_fit_pk('data/Isotope_dist_fit.pk')  # 1.寻找同位素峰，并返回记录不同价态的字典列表
        t0 = time.time()
        filter_MZ, max_int= pre_process_data(filter_MZ, max_int,fit_list, ppm, bound_th,
                                                     bound_intensity)
        # save_filter_MZ(filter_MZ, filter_MZ_path)  # 保存过滤结果
        t1 = time.time()
        print(t1 - t0)
    return filter_MZ, max_int


if __name__ == '__main__':
    filter_MZ_path = "data/filter_MZ.mgf"
    spectra_path = "data/938.mgf"
    candidate_path = "data/HP_20_mass.csv"
    # 单位ppm,ppm常用10
    candidates = input_candidate(candidate_path)
    filter_MZ, max_int = get_filter_MZ_pk(filter_MZ_path, spectra_path, candidates, ppm=10)
