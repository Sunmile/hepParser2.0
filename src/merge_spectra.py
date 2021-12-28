import numpy as np
import pandas as pd
import pickle as pk
from src.theory_sp import get_the_sp
from src.enumerate import get_comp_pk
from time import time


def get_the_mz_index(min_mz,max_mz, the_spectra, digit=2):
    the_iso = the_spectra[0]
    n_component = len(the_iso)
    id = 0
    dict_id_mz = {}
    dict_mz_id = {}
    for i in range(n_component):
        one_comp = the_iso[i]
        for j in range(len(one_comp)):
            iso = one_comp[j]
            iso_mz = iso[0]
            for mz in iso_mz:
                if mz > max_mz or mz < min_mz:
                    continue
                dict_id_mz[id] = mz
                new_mz = mz * np.power(10, digit)
                key_mz = int(np.round(new_mz))
                if key_mz in dict_mz_id.keys():
                    dict_mz_id[key_mz].append(id)
                else:
                    dict_mz_id[key_mz] = [id]
                if (key_mz + 1) in dict_mz_id.keys():
                    dict_mz_id[key_mz + 1].append(id)
                else:
                    dict_mz_id[key_mz + 1] = [id]
                if (key_mz - 1) in dict_mz_id.keys():
                    dict_mz_id[key_mz - 1].append(id)
                else:
                    dict_mz_id[key_mz - 1] = [id]
                id += 1
    return dict_mz_id, dict_id_mz


def match_index(mz_list, int_list, dict_mz_id, dict_id_mz, digit=2, ppm=10):
    ncol = len(dict_id_mz.keys())
    nrow = len(mz_list)
    res_list = np.zeros([nrow,ncol],dtype=int)
    all_dict_list = []
    count=0
    for i in range(len(mz_list)):
        one_mz_list = mz_list[i]
        one_int_list = int_list[i]
        dict_i_test = {}
        for j in range(len(one_mz_list)):
            one_mz = one_mz_list[j]
            one_int = one_int_list[j]
            if one_int <=0:
                continue
            new_mz = np.round(one_mz*np.power(10,digit))
            if new_mz in dict_mz_id.keys():
                candi_the_id = dict_mz_id[new_mz]
                for cid in candi_the_id:
                    the_mz = dict_id_mz[cid]
                    if np.abs(the_mz-one_mz)*1E6/the_mz < ppm:
                        res_list[i][cid] = one_int
                        if cid in dict_i_test.keys():
                            dict_i_test[cid].append(one_int)
                            count+=1
                            print('count:',count,', Has found conflict peak:', cid, one_int)
                        else:
                            dict_i_test[cid] = [one_int]
        all_dict_list.append(dict_i_test)
    return res_list, all_dict_list


def filter_mz_list(res_list,th=0.7):
    # idx = np.argwhere(np.all(res_list[..., :] == 0, axis=0))
    # new_res_list = np.delete(res_list,idx,axis=1)
    df_res = pd.DataFrame(res_list)
    df_res.drop(columns=df_res.columns[((df_res==0).mean()>th)],axis=1,inplace=True)
    return df_res


def get_merged_sp(mz_list, int_list, the_spectra, th=0.7,ppm=10, digit=2):
    min_list = [np.min(x) for x in mz_list]
    max_list = [np.max(x) for x in mz_list]
    min_mz = np.min(min_list)
    max_mz = np.max(max_list)
    t1 = time()
    dict_mz_id, dict_id_mz = get_the_mz_index(min_mz,max_mz,the_spectra,digit)
    res_list, all_dict_list = match_index(mz_list, int_list,dict_mz_id,dict_id_mz,digit,ppm)
    t2 = time()
    print(res_list.shape)
    print('match time:',t2-t1)
    df_res= filter_mz_list(res_list,th)
    sel_index = df_res.columns
    for one_index in sel_index:
        for i in range(len(all_dict_list)):
            if one_index not in all_dict_list[i].keys():
                tmp = [0]
            else:
                tmp = all_dict_list[i][one_index]
            tmp = ','.join([str(x) for x in tmp])
            df_res.loc[i,one_index] = tmp
    t3 = time()
    ind_list = df_res.columns
    sel_mz_list = []
    for ind in ind_list:
        mz = dict_id_mz[ind]
        sel_mz_list.append(str(np.round(mz,4)))
    df_res.columns=sel_mz_list
    print('filter time:',t3-t2)
    # df_res.to_csv('./../data/save_merged_BEH_peak39_ppm50_'+str(th)+'test2.csv')
    print(df_res.shape)
    print('save time:',time()-t3)


if __name__ == '__main__':
    dir = './../data/BEH_test_merge_all.pk'
    # # dir = './../data/test_merge_12.pk'
    # min_dp = 2
    # max_dp = 2
    # ppm = 10
    # max_ion = 10
    # digit = 2
    # th = 0.8
    # is_HNa = 'HNa'
    # t1 = time()
    with open(dir, 'rb') as f:
        data = pk.load(f)
    t2 = time()
    # print('read data:',t2-t1)
    mz_list = data[0]
    int_list = data[1]
    # the_HP = get_comp_pk('./../data/enox_', min_dp, max_dp)
    # the_SP_path = './../data/enox_' + str(min_dp) + '_' + str(max_dp) + '_sp_0.1_noloss_'
    # the_spectra = get_the_sp(the_SP_path, the_HP[3], the_HP[2], 1, max_ion, is_HNa=is_HNa)

    # get_merged_sp(mz_list,int_list,the_spectra,th,ppm,digit)
    # print('total time:',time()-t1)
    ####################################################
    # count_com=[]
    # max_ion = 10
    # is_HNa = 'HNa'
    # for i in range(2,21):
    #     min_dp=i
    #     max_dp=i
    #     count1 = 0
    #     count2 = 0
    #     the_HP = get_comp_pk('./../data/enox_', min_dp, max_dp)
    #     the_SP_path = './../data/enox_' + str(min_dp) + '_' + str(max_dp) + '_sp_0.1_noloss_'
    #     the_spectra = get_the_sp(the_SP_path, the_HP[3], the_HP[2], 1, max_ion, is_HNa=is_HNa)
    #     count1 = len(the_HP[3])
    #     for j in range(count1):
    #         count2 +=len(the_spectra[0][j])
    #     count_com.append([i,count1,count2])
    #     print(count_com[-1])
    # df_count = pd.DataFrame(count_com)
    # print(df_count)
    # df_count.columns=['dp','comp','isotope']
    # df_count.to_csv('./../data/count_dp.csv',index=None)
    tmp_list = []
    for i in range(len(mz_list)):
        tmp_list.extend(mz_list[i])
    tmp_arr = np.array(tmp_list)
    tmp_arr_uni = np.unique(tmp_arr)
    print(len(tmp_arr_uni))