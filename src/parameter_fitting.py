import numpy as np
import pandas as pd
from sklearn.metrics import r2_score
from sklearn.linear_model import Lasso
from src.theory_sp import *
from src.enumerate import get_comp_pk


# get two dict
# dict_pos_site : key: (3,2), value: the site in list
# dict_site_pos: key: the site in list, value:(3,2)
def get_site_dict(n_iso_list):
    dict_pos_site={}
    dict_site_pos={}
    for i in range(len(n_iso_list)):
        for j in range(n_iso_list[i]):
            pos = (i,j)
            site = sum(n_iso_list[0:i])+j
            dict_pos_site[pos] = site
            dict_site_pos[site] = pos
    return dict_pos_site,dict_site_pos


def transform_sp(the_SP_path, the_HP, max_ion, is_HNa):
    dict_mass_int = {}
    mass_list = []
    the_spectra = get_the_sp(the_SP_path, the_HP[3], the_HP[2], 1, max_ion, is_HNa=is_HNa)
    the_iso = the_spectra[0]
    n_component = len(the_iso)
    n_iso_list = []
    for i in range(n_component):
        n_iso_list.append(len(the_iso[i]))
    n_all_iso = sum(n_iso_list)
    dict_pos_site, dict_site_pos = get_site_dict(n_iso_list)
    for i in range(n_component):
        n_iso = len(the_iso[i])
        for j in range(n_iso):
            iso_mz_int = the_iso[i][j]
            iso_mz = iso_mz_int[0]
            iso_int = iso_mz_int[1]
            pos = (i, j)
            site = dict_pos_site[pos]
            for k in range(len(iso_mz)):
                tmp_mz = iso_mz[k]
                tmp_int = iso_int[k]
                if tmp_mz not in dict_mass_int.keys():
                    mass_list.append(tmp_mz)
                    row = np.zeros(n_all_iso)
                else:
                    row = dict_mass_int[tmp_mz]
                row[site] = tmp_int
                dict_mass_int[tmp_mz] = row
    mass_list = sorted(mass_list)
    return dict_mass_int, mass_list, dict_pos_site, dict_site_pos


def get_original_feature(mass_list, dict_mass_int):
    fea = []
    for x in mass_list:
        row = dict_mass_int[x]
        fea.append(row)
    return fea


def get_combine_feature(mass_list, dict_mass_int, ppm=20):
    fea = []
    new_mass_list = []
    pre_mass = 0
    n=1
    for i in range(len(mass_list)):
        tmp_mz = mass_list[i]
        win = tmp_mz * ppm * 1E-6
        tmp_row = dict_mass_int[tmp_mz]
        if tmp_mz -pre_mass < win:
            new_mass = np.round((pre_mass*n +tmp_mz)/(n+1),4)
            new_row = fea[-1]+tmp_row
            n+=1
            new_mass_list[-1] = new_mass
            fea[-1] = new_row
            pre_mass = new_mass
        else:
            n=1
            new_mass_list.append(tmp_mz)
            fea.append(tmp_row)
            pre_mass = tmp_mz
    return fea, new_mass_list


def generate_label(mass_list, exp_sp, max_int, ppm=20,scale=100):
    label = []
    i = 0
    j = 0
    while i<len(mass_list) and j<len(exp_sp):
        tmp_mz = mass_list[i]
        exp_mz = exp_sp[j][0]
        exp_int = exp_sp[j][1]/max_int
        win = tmp_mz *ppm * 1E-6
        if abs(tmp_mz-exp_mz)< win:
            label.append(np.round(exp_int*scale,4))
            i +=1
            j +=1
        elif tmp_mz > exp_mz:
            j += 1
        elif tmp_mz < exp_mz:
            i += 1
            label.append(0)
    return label


def lasso_fit(feature, label, alpha=1):
    lasso = Lasso(alpha=alpha, positive=True, max_iter=5000)
    lasso.fit(feature, label)
    y_pred = lasso.predict(feature)
    r2_score_lasso = r2_score(label, y_pred)
    w = lasso.coef_
    print('R2 score:', r2_score_lasso)
    # if alpha == 0.5:
    #     pd.DataFrame(w).to_csv('result/w_0.5.csv', header=None, index=None)
    # pl(w)
    return w, r2_score_lasso


def get_fit_result(w, dict_site_pos):
    dict_pos_w={}
    dict_comp_w={}
    for i in range(len(w)):
        if w[i]>0:
            pos = dict_site_pos[i]
            tmp_comp = pos[0]
            dict_pos_w[pos] = w[i]
            if tmp_comp not in dict_comp_w.keys():
                dict_comp_w[tmp_comp] = w[i]
            else:
                dict_comp_w[tmp_comp] += w[i]
    return dict_pos_w,dict_comp_w


if __name__ == '__main__':
    exp_dir = './../data/F12_merge.csv'
    dp = 4
    the_HP = get_comp_pk('./../data/enox_',dp)
    the_SP_path = './../data/enox_' + str(dp) + '_sp_0.1_noloss_'
    dict_mass_int, mass_list, dict_pos_site, dict_site_pos = \
        transform_sp(the_SP_path, the_HP, max_ion=10, is_HNa=True)
    ori_fea = get_original_feature(mass_list,dict_mass_int)
    comb_fea,new_mass_list = get_combine_feature(mass_list,dict_mass_int)
    print('ori mass list len:',len(mass_list))
    print('comb mass list len:',len(new_mass_list))
    exp_sp = pd.read_csv(exp_dir)

    exp_sp = np.array(exp_sp)
    max_int = max(exp_sp[:, 1])
    label = generate_label(new_mass_list, exp_sp, max_int)
    w, r2_score_lasso = lasso_fit(comb_fea,label, alpha=0.01)
    dict_pos_w, dict_comp_w = get_fit_result(w, dict_site_pos)
    print('R2 score:', r2_score_lasso)
    for x in dict_comp_w.keys():
        comp = the_HP[2][x]
        tmp_w = dict_comp_w[x]
        print(comp, tmp_w)
















