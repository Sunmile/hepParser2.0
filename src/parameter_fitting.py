import numpy as np
import pandas as pd
from sklearn.metrics import r2_score
from sklearn.linear_model import Lasso
from src.theory_sp import get_the_sp
from src.enumerate import get_comp_pk
from src.utils import draw
from time import time

from scipy import sparse
import xlwt

from scipy.stats import poisson


# get two dict
# dict_pos_site : key: (3,2), value: the site in list
# dict_site_pos: key: the site in list, value:(3,2)
def get_site_dict(n_iso_list):
    dict_pos_site = {}
    dict_site_pos = {}
    for i in range(len(n_iso_list)):
        for j in range(n_iso_list[i]):
            pos = (i, j)
            site = sum(n_iso_list[0:i]) + j
            dict_pos_site[pos] = site
            dict_site_pos[site] = pos
    return dict_pos_site, dict_site_pos




# ## all isotopic peaks are independent
# def transform_sp(the_spectra):
#     dict_mass_int = {}
#     mass_list = []
#     the_iso = the_spectra[0]
#     n_component = len(the_iso)
#     n_iso_list = []
#     for i in range(n_component):
#         n_iso_list.append(len(the_iso[i]))
#     n_all_iso = sum(n_iso_list)
#     dict_pos_site, dict_site_pos = get_site_dict(n_iso_list)
#     for i in range(n_component):
#         n_iso = len(the_iso[i])
#         for j in range(n_iso):
#             iso_mz_int = the_iso[i][j]
#             iso_mz = iso_mz_int[0]
#             iso_int = iso_mz_int[1]
#             pos = (i, j)
#             site = dict_pos_site[pos]
#             for k in range(len(iso_mz)):
#                 tmp_mz = iso_mz[k]
#                 tmp_int = iso_int[k]
#                 if tmp_mz not in dict_mass_int.keys():
#                     mass_list.append(tmp_mz)
#                     row = np.zeros(n_all_iso,dtype=np.int8)
#                 else:
#                     row = dict_mass_int[tmp_mz]
#                 row[site] = tmp_int
#                 dict_mass_int[tmp_mz] = row
#     mass_list = sorted(mass_list)
#     return dict_mass_int, mass_list, dict_pos_site, dict_site_pos


# ## the isotopic peaks with same charge of one component are dependent
# def transform_sp(the_spectra, lamb):
#     dict_mass_int = {}
#     mass_list = []
#     the_iso = the_spectra[0]
#     charge_iso = the_spectra[2]
#     hna_iso = the_spectra[4]
#     n_component = len(the_iso)
#     n_iso_list = []
#     for i in range(n_component):
#         n_iso_list.append(len(np.unique(charge_iso[i])))
#     n_all_iso = sum(n_iso_list)
#     dict_pos_site, dict_site_pos = get_site_dict(n_iso_list)
#     for i in range(n_component):
#         n_iso = len(the_iso[i])
#         for j in range(n_iso):
#             tmp_charge = charge_iso[i][j]
#             iso_mz_int = the_iso[i][j]
#             iso_mz = iso_mz_int[0]
#             iso_int = iso_mz_int[1]
#             pos = (i, tmp_charge-1)
#             site = dict_pos_site[pos]
#             Na_n = hna_iso[i][j][1]
#             for k in range(len(iso_mz)):
#                 tmp_mz = iso_mz[k]
#                 p = poisson.pmf(Na_n,lamb)
#                 tmp_int = iso_int[k]*p
#                 if tmp_mz not in dict_mass_int.keys():
#                     mass_list.append(tmp_mz)
#                     row = np.zeros(n_all_iso,dtype=np.int8)
#                 else:
#                     row = dict_mass_int[tmp_mz]
#                 row[site] = tmp_int
#                 dict_mass_int[tmp_mz] = row
#     mass_list = sorted(mass_list)
#     return dict_mass_int, mass_list, dict_pos_site, dict_site_pos

#  transform sp to feature with filter
def transform_sp(the_spectra, the_HP, dict_exp, digit=2):
    dict_mass_int = {}
    # consider the filtered isotope would change the site in feature matrix, record the map of sites.
    dict_fea_site = {}
    col = 0
    mass_list = []
    the_iso = the_spectra[0]
    n_component = len(the_iso)
    n_iso_list = []
    for i in range(n_component):
        n_iso_list.append(len(the_iso[i]))
    n_all_iso = sum(n_iso_list)
    dict_pos_site, dict_site_pos = get_site_dict(n_iso_list)
    for i in range(n_component):
        comp = the_HP[2][i]
        if comp[0]==0 and comp[5]==1: # filter the invalid components
            continue
        n_iso = len(the_iso[i])
        for j in range(n_iso):
            iso_mz_int = the_iso[i][j]
            iso_mz = iso_mz_int[0]
            iso_int = iso_mz_int[1]
            pos = (i, j)
            site = dict_pos_site[pos]
            int_order = np.argsort(iso_int)
            tmp_iso_mz = np.array(iso_mz)[int_order]
            if not filter_iso(tmp_iso_mz[::-1], dict_exp, digit):
                continue
            dict_fea_site[col] = site
            for k in range(len(iso_mz)):
                tmp_mz = iso_mz[k]
                tmp_int = iso_int[k]
                if tmp_mz not in dict_mass_int.keys():
                    mass_list.append(tmp_mz)
                    row = np.zeros(n_all_iso,dtype=np.int8)
                else:
                    row = dict_mass_int[tmp_mz]
                row[col] = tmp_int
                dict_mass_int[tmp_mz] = row
            col += 1
    mass_list = sorted(mass_list)
    print('Has filtered iso count:', n_all_iso-col)
    return dict_mass_int, mass_list, dict_pos_site, dict_site_pos, dict_fea_site, col


def get_original_feature(mass_list, dict_mass_int, col):
    fea = []
    for x in mass_list:
        row = dict_mass_int[x]
        dict_mass_int[x] = row[:col]
        fea.append(row[:col])
    return fea


def get_combine_feature(mass_list, dict_mass_int, ppm=20):
    fea = []
    new_mass_list = []
    pre_mass = 0
    n = 1
    for i in range(len(mass_list)):
        tmp_mz = mass_list[i]
        win = tmp_mz * ppm * 1E-6
        tmp_row = dict_mass_int[tmp_mz]
        if tmp_mz - pre_mass < win:
            new_mass = np.round((pre_mass * n + tmp_mz) / (n + 1), 4)
            new_row = fea[-1] + tmp_row
            n += 1
            new_mass_list[-1] = new_mass
            fea[-1] = new_row
            pre_mass = new_mass
        else:
            n = 1
            new_mass_list.append(tmp_mz)
            fea.append(tmp_row)
            pre_mass = tmp_mz
    # fea = sparse.csr_matrix(fea)
    return fea, new_mass_list


def generate_label(mass_list, exp_sp, max_int, ppm=20, scale=100):
    label = []
    diff = []
    fit_mz = []
    i = 0
    j = 0
    tmp_diff = -1
    while i < len(mass_list) and j < len(exp_sp):
        tmp_mz = mass_list[i]
        # if np.abs(tmp_mz -510.66)<0.05:
        #     a=1
        exp_mz = exp_sp[j][0]
        exp_int = exp_sp[j][1] / max_int
        win = tmp_mz * ppm * 1E-6
        if abs(tmp_mz - exp_mz) < win:
            if tmp_diff == -1:
                label.append(np.round(exp_int * scale, 4))
                fit_mz.append(exp_mz)
                tmp_diff = abs(tmp_mz - exp_mz)
                j += 1
            elif tmp_diff > abs(tmp_mz - exp_mz):
                label[-1] = np.round(exp_int * scale, 4)
                fit_mz[-1] = exp_mz
                tmp_diff = abs(tmp_mz - exp_mz)
                j += 1
            else:
                i += 1
                diff.append(tmp_diff)
                tmp_diff = -1

        elif tmp_mz > exp_mz:
            j += 1
        elif tmp_mz < exp_mz:
            if tmp_diff == -1:
                label.append(0)
                fit_mz.append(tmp_mz)
            else:
                diff.append(tmp_diff)
                tmp_diff = -1
            i += 1

    if len(diff)< len(label):
        diff.append(tmp_diff)
    if len(label)<len(mass_list):
        label.extend([0]*(len(mass_list)-len(label)))
        fit_mz.extend(mass_list[len(label):-1])
    return label,fit_mz


def lasso_fit(feature, label, alpha=0.01):
    lasso = Lasso(alpha=alpha, positive=True, tol=1E-6, random_state=0, max_iter=50000)
    lasso.fit(feature, label)
    y_pred = lasso.predict(feature)
    r2_score_lasso = r2_score(label, y_pred)
    w = lasso.coef_
    # print('R2 score:', r2_score_lasso)
    # if alpha == 0.5:
    #     pd.DataFrame(w).to_csv('result/w_0.5.csv', header=None, index=None)
    # pl(w)
    return w, r2_score_lasso


def get_fit_result(w, dict_site_pos, dict_fea_site):
    dict_pos_w = {}
    dict_comp_w = {}
    for i in range(len(w)):
        if w[i] > 0:
            site = dict_fea_site[i]
            pos = dict_site_pos[site]
            tmp_comp = pos[0]
            dict_pos_w[pos] = w[i]
            if tmp_comp not in dict_comp_w.keys():
                dict_comp_w[tmp_comp] = w[i]
            else:
                dict_comp_w[tmp_comp] += w[i]
    return dict_pos_w, dict_comp_w


def get_exp_dict(exp_sp, digit=2):
    dict_exp = {}
    for i in range(len(exp_sp)):
        mz = exp_sp[i,0]*np.power(10,digit)
        m = int(np.round(mz))
        inten = exp_sp[i,1]
        if inten>0:
            dict_exp[m] = inten
    return dict_exp


# here only the iso which contain two continuous matched peaks would be keep
def filter_iso(iso_mz, dict_exp, digit=2):
    count = 0
    for i in range(len(iso_mz)):
        mz = iso_mz[i]*np.power(10,digit)
        m = int(np.round(mz))
        if m in dict_exp.keys() or m+1 in dict_exp.keys() or m-1 in dict_exp.keys():
            count += 1
        else:
            break
    if count >= 2:
        return True
    else:
        return False


def trans_show_format(comp):
    str_comp = '['
    if comp[0] == 1:
        str_comp='dA['
    for x in comp[1:5]:
        str_comp += str(x)+','
    str_comp = str_comp[:-1]
    str_comp += ']'
    if comp[5] ==1:
        str_comp += 'aG'
    if comp[6] ==1:
        str_comp += 'aM'
    return str_comp


def run_main(save_dir, exp_dir):
    # exp_dir = './../data/fitting/L8_orginal_peak.csv'
    lamb = 0.8  # 1.5 for dp4, 0.8 for dp6
    min_dp = 10
    max_dp = 12
    # ppm0 = 50
    alpha = 0.01
    ppm = 20
    max_ion = 10
    digit = 2
    is_HNa = 'HNa'
    exp_sp = pd.read_csv(exp_dir)
    exp_sp = np.array(exp_sp)
    max_int = max(exp_sp[:, 1])
    the_HP = get_comp_pk('./../data/enox_', min_dp, max_dp)
    the_SP_path = './../data/enox_' + str(min_dp) + '_' + str(max_dp) + '_sp_0.1_noloss_'
    the_spectra = get_the_sp(the_SP_path, the_HP[3], the_HP[2], 1, max_ion, is_HNa=is_HNa)
    dict_exp = get_exp_dict(exp_sp,digit=digit)
    dict_mass_int, mass_list, dict_pos_site, dict_site_pos, dict_fea_site,col =\
        transform_sp(the_spectra,the_HP,dict_exp,digit)
    # dict_mass_int, mass_list, dict_pos_site, dict_site_pos = transform_sp(the_spectra, lamb)
    ori_fea = get_original_feature(mass_list, dict_mass_int,col)
    comb_fea, new_mass_list = get_combine_feature(mass_list, dict_mass_int, ppm)
    print('ori mass list len:', len(ori_fea), len(ori_fea[0]))
    print('comb mass list len:', len(comb_fea), len(comb_fea[0]))
    label,fit_mz = generate_label(new_mass_list, exp_sp, max_int, ppm=ppm, scale=100)
    w, r2_score_lasso = lasso_fit(comb_fea, label, alpha=alpha)
    dict_pos_w, dict_comp_w = get_fit_result(w, dict_site_pos, dict_fea_site)
    print('R2 score:', r2_score_lasso)
    comp_w = []
    for x in dict_comp_w.keys():
        comp = the_HP[2][x]
        tmp_w = dict_comp_w[x]
        com_str = trans_show_format(comp)
        comp_w.append([com_str,np.round(tmp_w, 5)])
        print(com_str, '\t', np.round(tmp_w, 5))
    comp_result = pd.DataFrame(comp_w,columns=['Components','Weight'])
    comp_result.sort_values(by='Weight', ascending=False, inplace=True)
    annotate = []
    for i in range(len(new_mass_list)):
        exp_int = label[i]
        if exp_int > 0:
            mz = new_mass_list[i]
            exp_mz = fit_mz[i]
            the_int = np.dot(comb_fea[i], w)
            one_fea = np.array(comb_fea[i]) * np.array(w)
            fea_site = np.where(one_fea > 0)[0]
            poss = []
            HNas = []
            comps = []
            mzs = []
            ws = []
            ppms = []
            for x in fea_site:
                site = dict_fea_site[x]
                pos = dict_site_pos[site]
                HNa = the_spectra[4][pos[0]][pos[1]]
                comp = the_HP[2][pos[0]]
                com_str = trans_show_format(comp)
                tmp_mz = the_spectra[0][pos[0]][pos[1]][0][0]
                one_w = np.round(dict_pos_w[pos], 5)
                ppm = np.round(np.abs((tmp_mz-exp_mz))*1E6/tmp_mz,2)
                if ppm > 50:
                    ppm = -1
                mzs.append(tmp_mz)
                comps.append(com_str)
                poss.append(pos)
                HNas.append(HNa)
                ws.append(one_w)
                ppms.append(ppm)
            annotate.append([mz, exp_mz, exp_int, the_int, poss, comps, mzs, ws, HNas, ppms])
    annotate = pd.DataFrame(annotate, columns=['mz', 'exp_mz', 'exp_int', 'the_int', 'position', 'components', 'the_mz', 'weight', 'HNas','diff_ppm'])
    # annotate.to_csv(save_dir, index=None)
    writer = pd.ExcelWriter(save_dir)
    annotate.to_excel(writer, sheet_name='annotation', index=None)
    comp_result.to_excel(writer, sheet_name='components', index=None)
    writer.save()
    return r2_score_lasso, exp_sp


if __name__ == '__main__':
    # exp_dir = './../data/BEH_peak5.csv'
    start = time()
    # sample = 'X19'
    #
    # exp_dir = './../data/fitting/' + sample + '_merged.csv'
    # dir = './../data/fitting/' + sample + '_20_fit.xlsx'
    # fig_dir = './../data/figure/' + sample + '_20_fit.png'
    sample = 'Luna_RT 26.11-26.80_dp12'
    exp_dir = './../data/merge_dicp/'+sample+'.csv'
    dir ='./../data/merge_dicp/'+sample+'_10_fit.xlsx'
    fig_dir ='./../data/merge_dicp/' + sample + '_10_fit.png'

    r2_score_lasso, exp_sp = run_main(dir, exp_dir)

    # exp_dir = './../data/BEH_peak5.csv'
    # exp_sp = pd.read_csv(exp_dir)
    # exp_sp = np.array(exp_sp)
    # r2_score_lasso= 0.99973
    print('time: ', np.round(time()-start,2))
    draw(dir, r2_score_lasso, exp_sp, fig_dir)
    # test = [1,1,2,0,6,1,0]
    # print(trans_show_format(test))



