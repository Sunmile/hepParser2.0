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


def filter_spectra_with_same_int(MZ_int, th=100):
    dict_int = {}
    for i in range(len(MZ_int)):
        if MZ_int[i][1] in dict_int.keys():
            dict_int[MZ_int[i][1]] += 1
        else:
            dict_int[MZ_int[i][1]] = 1
    filter_MZ = []
    for i in range(len(MZ_int)):
        if MZ_int[i][1] == 0 and (len(filter_MZ) == 0 or filter_MZ[-1][1] != 0):
            filter_MZ.append(MZ_int[i])
            continue
        if dict_int[MZ_int[i][1]] <= th:
            filter_MZ.append(MZ_int[i])
    return filter_MZ


def format_mz(peaks):
    filter_MZ = []
    for i in range(len(peaks)):
        if peaks[i][0] >= 0:  # 350
            filter_MZ.append([np.round(peaks[i][0], 4), peaks[i][1]])
    return filter_MZ


def is_min(peaks,idx, win=5): ## 判断给定idx的峰是否是i，i+win中最小的
    end = len(peaks)-1
    if idx+win <end:
        end = idx+win
    min_peaks = peaks[idx][1]
    for i in range(idx+1, end+1):
        if peaks[i][1] > 3 * min_peaks:
            return True
        if min_peaks > peaks[i][1]:
            return False
    return True


def merge_peaks_simple(peaks, max_int=0, threshold=0.4):
    signal_list = []
    for i in range(0,len(peaks)-1):
        signal_list.append(peaks[i+1][1]-peaks[i][1])
    signal_list.append(0)
    new_peaks = []
    tmp_max = 0
    tmp_mz = []
    tmp_intensity = []
    for i in range(len(peaks)):
        if tmp_max<peaks[i][1]:
            tmp_max = peaks[i][1]
        if signal_list[i-1]<=0 and not signal_list[i]<=0:
            if is_min(peaks,i,3):
                tmp_int = peaks[i][1]
                if tmp_int <= tmp_max * threshold:
                    new_mz = np.sum(np.multiply(np.array(tmp_mz),np.array(tmp_intensity)))/np.sum(tmp_intensity)
                    new_int = np.sum(tmp_intensity)
                    if max_int < new_int:
                        max_int = new_int
                    if not np.isnan(new_mz):
                        new_peaks.append([np.round(new_mz,5),new_int])
                    tmp_mz = []
                    tmp_intensity = []
                    tmp_max = 0
                else:
                    tmp_mz.append(peaks[i][0])
                    tmp_intensity.append(peaks[i][1])
            else:
                tmp_mz.append(peaks[i][0])
                tmp_intensity.append(peaks[i][1])
        else:
            tmp_mz.append(peaks[i][0])
            tmp_intensity.append(peaks[i][1])
    for i in range(len(new_peaks)):
        new_peaks[i].append(new_peaks[i][1] / max_int)
    return new_peaks, max_int


## 新的过滤方法 测试
def get_mean_of_win(sp_df, win=1, step=0.5):
    sp_df['key1'] = sp_df['mz']//win
    sp_df['key2'] = (sp_df['mz']+0.5)//win -step
    dict_mean = {}
    for x in sp_df['key1'].unique():
        m = np.mean(sp_df[sp_df['key1']==x]['Int'])
        dict_mean[x]=m
    for x in sp_df['key2'].unique():
        m = np.mean(sp_df[sp_df['key2']==x]['Int'])
        dict_mean[x]=m
    return sp_df,dict_mean


def find_min_mean(dict_mean, start, win=5):
    min_mean = -1
    for i in np.arange(start-win/2, start+win/2+0.1, 0.5):
        if i in dict_mean.keys():
            tmp = dict_mean[i]
            if min_mean < 0 or tmp < min_mean:
                min_mean = tmp
    if min_mean < 0:
        min_mean = 0
    return min_mean


def label_noise(sp_df, dict_mean, win=5):
    noise_th = []
    sp_df['noise'] = 0
    for x in sp_df['key1'].unique():
        # print('x',x)
        th = find_min_mean(dict_mean,x)
        noise_th.append([x, th])
        sp_df.loc[sp_df[sp_df['key1']==x][sp_df['Int']<=th].index.tolist(),'noise']=1
    return sp_df, np.array(noise_th)


# 预处理的整个流程：峰合并、峰过滤、校准谱
def pre_process_data(filter_MZ, max_int, fit_list, ppm=10, bound_th=0.001, bound_intensity=300):
    '''
    # baseline
    filter_MZ = filter_spectra(filter_MZ, 0.01, max_int, 100)  # 过滤谱文件，去除强度小于100且相对丰度低于0.01的峰
    filter_MZ = local_denoise(filter_MZ, 0.1, 0.2)  # 针对局部进行峰过滤，只保留局部最高峰
    num_matched = match_peak_ppm(filter_MZ, max_int, candidates, fit_list, ppm)
    print("baseline:", "max_int=" + str(max_int), "peak_num=" + str(len(filter_MZ)), "matched_num=" + str(num_matched))
    '''
    th = 1000
    filter_MZ = format_mz(filter_MZ)
    save_test = pd.DataFrame(filter_MZ)
    save_test.to_csv('data/orginal_peak.csv', header=True, index=None)
    # 原始峰
    # num_matched = match_peak_ppm(filter_MZ, max_int, candidates, fit_list, ppm)
    # print("origin:", "peak_num=" + str(len(filter_MZ)), "matched_num=" + str(num_matched))
    # 初始过滤机器噪音
    filter_MZ = filter_spectra_with_same_int(filter_MZ, th)
    save_test = pd.DataFrame(filter_MZ)
    save_test.to_csv('data/filter1_peak.csv', header=True, index=None)
    # 峰合并
    # filter_MZ, max_int = merge_peaks_with_original_mz(filter_MZ, max_int, 0.6)  # 2.针对局部进行峰过滤，只保留局部最高峰
    # filter_MZ, max_int = merge_peaks(filter_MZ, max_int, 0.2)  # 2.针对局部进行峰过滤，只保留局部最高峰
    filter_MZ, max_int = merge_peaks_simple(filter_MZ, max_int, 0.4)
    save_test = pd.DataFrame(filter_MZ)
    save_test.to_csv('data/merged.csv', header=True, index=None)
    # num_matched = match_peak_ppm(filter_MZ, max_int, candidates, fit_list, ppm)
    # print("merge:", "max_int=" + str(max_int), "peak_num=" + str(len(filter_MZ)), "matched_num=" + str(num_matched))

    # 过滤

    df_B = pd.DataFrame(filter_MZ)
    df_B.columns = ['mz', 'Int', 'Int_R']
    df_B, dict_mean = get_mean_of_win(df_B)
    df_B, noise_th = label_noise(df_B, dict_mean)
    df_noise = df_B[df_B['noise'] == 1]
    df_denoise = df_B[df_B['noise']==0]
    filter_MZ = np.array(df_denoise[['mz','Int','Int_R']]).tolist()

    # filter_MZ = filter_spectra(filter_MZ, bound_th, max_int, bound_intensity)  # 3.过滤谱文件，去除强度小于300且相对丰度低于0.001的峰
    save_test = pd.DataFrame(filter_MZ)
    save_test.to_csv('data/final_filter.csv', header=True, index=None)
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
        fit_list = get_fit_pk('../data/Isotope_dist_fit.pk')  # 1.寻找同位素峰，并返回记录不同价态的字典列表
        t0 = time.time()
        filter_MZ, max_int = pre_process_data(filter_MZ, max_int, fit_list, ppm, bound_th,
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
