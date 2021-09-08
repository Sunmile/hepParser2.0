from src.peak_align import *
from src.isotopic_detect import *
from src.match_score import *
from src.label import *
from src.significance_test import *

def get_filter_MZ(origin_MZ, max_int, fit_list, the_HP, the_sp_path, max_ion=7, ppm=20, bound_th=0.001, bound_intensity=0, is_HNa='H'):
    from src.preprocess import pre_process_data
    # 因为2个文件相互引用了，所以在这里import

    s = time()
    top_n = 5
    isp_thresh = 0.5
    the_SP_path = the_sp_path
    t0 = time()
    print('read:', t0 - s)
    filter_MZ, max_int = pre_process_data(origin_MZ, max_int, fit_list, ppm, bound_th,
                                          bound_intensity)
    total_int = get_total_intensity(filter_MZ, max_int)
    # save_file(filter_MZ, 'data/filter_MZ')
    t1 = time()
    print('preprocess', t1 - t0)
    dict_score_list, dict_list, dict_the_iso_list = get_isotopic(filter_MZ, fit_list, ppm, max_ion, 5, isp_thresh)
    tmp_save = [dict_list, dict_the_iso_list]
    # save_pk(tmp_save,'./data/isotopic_draw.pk')

    t2 = time()
    print('find isotopic:', t2 - t1)

    the_spectra = get_the_sp(the_SP_path, the_HP[3],the_HP[2], 1, max_ion, is_HNa=is_HNa)
    t7 = time()
    print('generate the sp:', t7 - t2)
    dict_list, delta = align_peak(dict_list, the_HP, the_spectra, top_n, ppm)
    exp_isp = get_exp_isp(dict_list, max_int)
    peaks = filter_MZ
    for i in range(0, len(peaks)):
        peaks[i][0] = np.round(peaks[i][0] + delta, 5)
    filter_MZ = peaks
    return filter_MZ, max_int, total_int, exp_isp, the_spectra, dict_list


def data_process(peaks, exp_isp, max_int, total_int, the_spectra, dict_list, the_HP, chb_dA, chb_aM, chb_aG, min_dp=0,max_dp=20, ppm=20):
    prob = 0.95 # 衍生峰得分权重系数
    confidence_list = get_probability(peaks, the_HP,the_spectra,ppm)
    all_comp, all_atom, all_comp_mass, all_comp_score, dict_exp_match, \
    dict_theo_list, all_score_list, all_score_mass, all_comp_lost, all_comp_z, all_comp_HNa = \
        calculate_all_hp_score(the_spectra, dict_list, the_HP, max_int, chb_dA=chb_dA, chb_aM=chb_aM, chb_aG=chb_aG,
                               min_dp=min_dp, max_dp=max_dp, ppm=ppm, prob=prob)
    t3 = time()
    # print('match time:', t3 - t7)
    match_result = [all_comp, all_atom, all_comp_mass, all_comp_score, dict_exp_match, dict_theo_list, all_score_list,
                    all_score_mass, all_comp_lost, all_comp_z, all_comp_HNa]
    # save_pk(match_result, result_path + 'one_match_result0.5.pk')
    result2_score, candidate_sort_index = get_most_power_candidate_score(match_result)
    t4 = time()
    print('pick:', t4 - t3)
    label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order = get_n_label(match_result,
                                                                                          total_int,
                                                                                          candidate_sort_index,
                                                                                          len(candidate_sort_index))
    right_comp = []
    matched_count = []
    for i in range(len(key_with_order)):
        tmp_key = key_with_order[i]
        comp = str(dict_mass_comp[tmp_key])
        tmp_c = len(dict_mass_family[tmp_key])
        matched_count.append(tmp_c)
        right_comp.append(comp)
    df_conf = pd.DataFrame(confidence_list, columns=['comp','prob', 'peak_count'])
    # df_conf.to_csv('result/df_conf.csv',index=None)
    df_conf['comp'] = [str(x) for x in df_conf['comp']]
    df_conf = df_conf[df_conf.comp.isin(right_comp)]
    if len(right_comp)<1:
        df_conf['log_p'] = []
    else:
        df_conf = get_pvalue(right_comp,matched_count,df_conf)
    # label_save = pd.DataFrame(np.array(label))
    # label_save.to_csv('result/label_save.csv', index=None)
    t5 = time()
    print('label:', t5 - t4)

    e = time()
    # print('Time cost', e - s)
    # label_info = get_n_result(match_result, max_candi + 1)

    label_info = [label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order, df_conf]
    sup_n = 200
    if len(result2_score[1]) > sup_n:
        label_info = [label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order[:sup_n], df_conf]
        return [match_result, label_info, sup_n, result2_score[1][:sup_n]]
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
    from src.preprocess import get_filter_MZ_pk, pre_process_data

    s = time()
    is_HNa=True
    ppm = 20
    top_n = 5
    prob = 0.95
    isp_thresh = 0.8
    bound_th = 0.001
    bound_intensity = 0
    the_HP = get_comp_pk('data/delta',6)  # 枚举的理论肝素结构
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
    dict_score_list, dict_list, dict_the_iso_list = get_isotopic(filter_MZ, fit_list, ppm, 7, 5,
                                                                    isp_thresh)  # 寻找同位素峰，并返回记录不同价态的字典列表
    t2 = time()
    print('find isotopic:', t2 - t1)
    the_spectra = get_the_sp(the_SP_path, the_HP[3],the_HP[2], 1,max_ion=7,is_HNa=is_HNa)
    t7 = time()
    print('generate the sp:', t7 - t2)
    dict_list, delta = align_peak(dict_list, the_HP, the_spectra, top_n, ppm)
    all_comp, all_atom, all_comp_mass, all_comp_score, dict_exp_match, \
    dict_theo_list, all_score_list, all_score_mass, all_comp_lost, all_comp_z, all_comp_HNa = \
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
