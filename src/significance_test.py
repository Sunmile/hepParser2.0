import numpy as np
import pandas as pd
from scipy.stats import binom
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector


dict_binom = {}


## 返回每个成分随机被匹配上一个峰的概率
def get_probability(peaks, the_HP, the_spectra,ppm=20):
    min_mz = np.min(np.array(peaks)[:,0])
    max_mz = np.max(np.array(peaks)[:,0])
    all_len = max_mz-min_mz
    confidence_list = []
    for i in range(len(the_HP[2])):
        comp = the_HP[2][i]
        iso_list = the_spectra[0][i]
        iso_len = 0
        count = 0
        for j in range(len(iso_list)):
            first_mz = iso_list[j][0][0]
            if first_mz > min_mz and first_mz < max_mz:
                win = first_mz * ppm * 0.000001
                iso_len += win
                count += 1
        p = iso_len/all_len
        confidence_list.append([comp, p, count])
    confidence_list = np.array(confidence_list)
    # R_stats = importr('stats')
    # p_adj = R_stats.p_adjust(FloatVector(confidence_list[:,1]), method='BY')
    # p_adj2 = R_stats.p_adjust(FloatVector(confidence_list[:,1]), method='fdr')
    # confidence_list = np.hstack((confidence_list,np.array(p_adj).reshape(-1,1)))
    # confidence_list = np.hstack((confidence_list, np.array(p_adj2).reshape(-1, 1)))
    # df = pd.DataFrame(confidence_list)
    # df.to_csv('result/confidence.csv', index=None)
    return confidence_list



def get_binom_pmf(k,n,p):
    if (k,n,p) in dict_binom.keys():
        return dict_binom[(k,n,p)]
    if (k,n-1,p) in dict_binom.keys():
        tmp= dict_binom[(k,n-1,p)]
        res = tmp*(1-p)*n/(n-k)
        dict_binom[(k,n,p)]=res
        return res
    if (k-1,n,p) in dict_binom.keys():
        tmp=dict_binom[(k-1,n,p)]
        res = tmp*p/(1-p)*(n-k+1)/k
        dict_binom[(k,n,p)]=res
        return res
    res = binom.pmf(k,n,p)
    dict_binom[(k,n,p)]=res
    return res


def get_p_sig(k, n, p):
    tmp_p = 0
    for i in range(k,n):
        tmp_p += get_binom_pmf(i,n,p)
    return tmp_p


def get_pvalue(right_comp, matched_count, df_conf):
    df_match = pd.DataFrame(np.array([right_comp,matched_count]).T,columns=['comp','matched'])
    df_conf = df_conf.merge(df_match,on='comp')
    pvalue = []
    df_conf.index = range(df_conf.shape[0])
    for i in range(df_conf.shape[0]):
        comp = df_conf.loc[i,'comp']
        prob = df_conf.loc[i,'prob']
        n = np.int(df_conf.loc[i,'peak_count'])
        k = np.int(df_conf.loc[i,'matched'])
        one_p = get_p_sig(k,n,prob)
        pvalue.append(one_p)
    df_conf['pvalue'] = pvalue
    R_stats = importr('stats')
    p_adj = R_stats.p_adjust(FloatVector(df_conf['pvalue']), method='BY')
    p_adj2 = R_stats.p_adjust(FloatVector(df_conf['pvalue']), method='fdr')
    df_conf['p1']=p_adj
    df_conf['p2']=p_adj2
    df_conf['log_p'] = -np.log10(df_conf['p2'])
    # df_conf.to_csv('result/conf_fine.csv',index=None)
    return df_conf
