import numpy as np
from src.utils import *
from src.theory_sp import *
from time import time


def score_match_isp(isp_1, isp_2, ppm):
    '''
    input:
    isp_1: 同位素峰簇1   理论
    isp_2: 同位素峰簇2   实际
    example: [[m/z1,m/z2,m/z3,m/z4,m/z5],[int1,int2,int3,int4,int5]]
    注意: isp_1, isp_2 的m/z 值必须对应相等
    output:

    score: 匹配得分
    '''
    mz_list1 = isp_1[0]
    mz_list2 = isp_2[0]
    int_list1 = isp_1[1]
    int_list2 = isp_2[1]
    int_weight = 0
    for x in int_list2:
        int_weight += np.log(x + 1)
        # int_weight += x
    # int_weight = np.log(np.sum(int_list2)+1)
    score_mz = 0
    for i in range(len(mz_list1)):
        diff = abs(mz_list1[i] - mz_list2[i])
        win = ppm * mz_list1[i] * 0.000001
        if diff <= win:
            score_mz += diff / win + 0.1
    # score_int = cos_sim(int_list1,int_list2)
    score_int = int_weight * (1 - get_JS(int_list1, int_list2, 0.5, 0.5))
    # score_int = get_KL(int_list1,int_list2)
    # score_all = score_mz * score_int
    score_all = score_int
    return score_all


def match_isotopic(exp_isp, the_isp, ppm, thresh, min_match_num, the_lost_list, the_z_list, the_HNa_list, p=0.9):
    '''
    input:
    exp_isp: 实验谱，形式是同位素峰簇的集合，每个同位素峰簇包含五个峰的m/z值，以及对应的int值。三重list，
    example: [[[m/z1,m/z2,m/z3,m/z4,m/z5],[int1,int2,int3,int4,int5]],[[m/z1,m/z2,m/z3,m/z4,m/z5],[int1,int2,int3,int4,int5]]]
    the_isp: 理论谱，格式与exp_isp 相同
    ppm: win = mass*ppm*10^-6
    win: 匹配窗口的大小，两个峰m/z 差值在win之内，被认为是匹配成功
    thresh: 单个同位素峰匹配得分阈值
    min_match_num: 	单个同位素峰簇中至少匹配的峰数


    output:
    global_score: 匹配的全局得分, 单次匹配的加和
    matched_isp_list: 被匹配上的同位素峰簇的集合
    ret_the_isp_list: 匹配上的同位素峰簇的对应理论峰簇的集合
    score_list: 单个同位素峰簇得分list
    score_mass: 匹配的单同位素峰

    '''

    '''
    执行操作:
    对于exp_isp 中的每一个 同位素峰簇，分别与 the_isp 中的同位素峰簇匹配
    若有不低于min_match_num个峰被匹配上，则认为该两个同位素峰簇匹配上了。
    匹配上的话，以 the_isp中的同位素峰簇为准，将exp_isp中的同位素峰簇补全，补充的峰的int值为0(注意，这步不能修改原exp_isp值),
    然后调用score_match_isp(isp_1,isp2,win) 获得 匹配得分。
    example1：
    exp_isp中：isp_1=[[1, 1.5, 2, 2.5, 3],[10, 20, 11, 10, 5]]
    the_isp中：isp_2=[[1, 2, 3, 4, 5],[20, 5, 10, 22, 4]]
    则 isp_1 更改为：[[1, 2, 3, 4, 5],[10, 11, 5, 0, 0]]
    example2：
    exp_isp中：isp_1=[[3, 4, 5, 6, 7],[10, 20, 11, 10, 5]]
    the_isp中：isp_2=[[1, 2, 3, 4, 5],[20, 5, 10, 22, 4]]
    则 isp_1 更改为：[[1, 2, 3, 4, 5],[0, 0, 10, 20, 11]]

    只有匹配得分超过thresh的，才会被加到matched_isp_list中, 加的是更改后的exp_isp 中的同位素峰簇。
    '''
    # 找到最长的间距
    MAX_L = 0
    for item in exp_isp:
        temp = item[0][-1] - item[0][0]
        if temp > MAX_L:
            MAX_L = temp

    for item in the_isp:
        temp = item[0][-1] - item[0][0]
        if temp > MAX_L:
            MAX_L = temp

    matched_isp_list = []
    matched_score_list = []
    global_score = 0
    # 记录matched_isp_list中每组匹配的实验谱和理论普索引，用于后期筛选
    index_exp = []
    index_the = []
    i = 0
    j0 = 0
    while (i < len(exp_isp)):
        # 对于每个实验谱中的峰簇，找到匹配分数最大的理论谱
        exp_item = exp_isp[i]
        while (exp_item[0][0] > the_isp[j0][0][0] + MAX_L):
            j0 = j0 + 1
            if j0 == len(the_isp):
                j0 = j0 - 1
                break

        if j0 == len(the_isp) - 1 and exp_item[0][0] > the_isp[j0][0][0] + MAX_L:
            break

        max_score = 0
        max_list = None
        max_j = j0
        j = j0
        while (j < len(the_isp)):
            the_item = the_isp[j]
            tmp_lost_sum = np.sum(the_lost_list[j])
            if the_item[0][0] > exp_item[0][0] + MAX_L:
                break

            counter = 0
            # 计算两个峰簇的匹配数量,并整理出补全后的峰簇temp_exp
            temp_exp = []
            temp_exp.append([])
            temp_exp.append([])
            match_str = ''
            win = ppm * the_item[0][0] * 0.000001
            for aim_mz in the_item[0]:
                # 两个峰簇匹配时，使用理论谱峰簇的m/z值在实验谱中搜索，如果找到win之内的，则记录，否则记0
                temp_exp[0].append(aim_mz)
                temp_value = 0
                k = 0
                match_flag = False
                for k in range(counter, len(exp_item[0])):
                    exp_mz = exp_item[0][k]
                    if abs(exp_mz - aim_mz) <= win and exp_item[1][k] != 0:
                        temp_value = exp_item[1][k]
                        temp_exp[0][-1] = exp_item[0][k]
                        counter = counter + 1
                        match_str += '1'
                        match_flag = True
                        break
                temp_exp[1].append(temp_value)
                if not match_flag:
                    match_str += '0'
            # 匹配数目足够，则将该匹配看做备选项
            if counter < min_match_num or match_str == '10101':
                j = j + 1
                continue

            # 计算备选项匹配得分
            this_score = score_match_isp(the_item, temp_exp, ppm)
            if tmp_lost_sum > 0:
                this_score = p * this_score

            if this_score <= thresh:
                j = j + 1
                continue

            # 备选项的分数比当前最大值大的话，更新最大值及相关索引
            if this_score > max_score:
                max_score = this_score
                max_list = temp_exp.copy()
                max_j = j
            j = j + 1

        # 若该实验峰簇未找到匹配项
        if max_score == 0:
            i = i + 1
            continue

        # 记录补全后的峰簇、累加匹配分数、记录两个索引
        matched_isp_list.append(max_list)
        matched_score_list.append(max_score)
        global_score = global_score + max_score
        index_exp.append(i)
        index_the.append(max_j)
        i = i + 1

    # 判断是否有重复
    if len(set(index_the)) == len(index_the):
        ret_the_isp_list = []
        ret_the_lost_list = []
        ret_the_z_list = []
        ret_the_HNa_list = []
        for i in index_the:
            ret_the_isp_list.append(the_isp[i])
            ret_the_lost_list.append(the_lost_list[i])
            ret_the_z_list.append(the_z_list[i])
            ret_the_HNa_list.append(the_HNa_list[i])
        score_mass = []  # 标记 得分的单同位素峰
        for i in range(len(matched_isp_list)):
            mass = matched_isp_list[i][0]
            score_mass.append(mass)
        return global_score, matched_isp_list, ret_the_isp_list, matched_score_list, score_mass, ret_the_lost_list, \
               ret_the_z_list, ret_the_HNa_list

    # 有重复的话
    index = []
    # 找到重复的项，存入做为set的index中
    for i in range(len(index_the)):
        j = i + 1
        while (j < len(index_the)):
            if (index_the[i] == index_the[j]):
                index_flag = 0
                for item in index:
                    if i in item:
                        item.add(j)
                        index_flag = 1
                        break
                    elif j in item:
                        item.add(i)
                        index_flag = 1
                        break
                if index_flag == 0:
                    temp_set = set()
                    temp_set.add(i)
                    temp_set.add(j)
                    index.append(temp_set)
            j = j + 1

    delete_counter = 0
    # 对于每一组重复
    for item in index:
        item_list = list(item)
        if len(item) == 2:
            i = item_list[0]
            j = item_list[1]
            # flag表征两类情况，1代表可合并，0代表不可合并应取分数最大
            flag = 0
            # if和elif代表两种可能出现合并的起始情况，若可合并，则最后flag会被置为1
            # if：i对应峰簇前半部分为0，j对应峰簇后半部分为0
            if matched_isp_list[i][1][0] == 0 and matched_isp_list[j][1][0] != 0:
                # 核心思路是找到i由0变非0，以及j由非0变0的位置，使用change做标记
                # 若上述改变只有一次，且未发生其他一些特殊情况，则认为i和j能够合并
                # elif里思路相同
                change = 0
                for k in range(len(matched_isp_list[i][1])):
                    if matched_isp_list[i][1][k] != 0 and matched_isp_list[j][1][k] != 0:
                        break
                    if change == 0 and matched_isp_list[i][1][k] == 0 and matched_isp_list[j][1][k] != 0:
                        continue
                    elif change == 0 and matched_isp_list[j][1][k] == 0 and matched_isp_list[i][1][k] == 0:
                        if k == len(matched_isp_list[i][1]) - 1:
                            flag = 1
                        continue
                    elif change == 0 and matched_isp_list[j][1][k] == 0 and matched_isp_list[i][1][k] != 0:
                        change = 1
                        if k == len(matched_isp_list[i][1]) - 1:
                            flag = 1
                        continue
                    elif change == 1 and matched_isp_list[i][1][k] == 0 and matched_isp_list[j][1][k] == 0:
                        if k == len(matched_isp_list[i][1]) - 1:
                            flag = 1
                        continue
                    elif change == 1 and matched_isp_list[i][1][k] != 0 and matched_isp_list[j][1][k] == 0:
                        if k == len(matched_isp_list[i][1]) - 1:
                            flag = 1
                        continue
                    elif change == 1 and matched_isp_list[i][1][k] == 0 and matched_isp_list[j][1][k] != 0:
                        break

            # elif：i对应峰簇后半部分为0，j对应峰簇前半部分为0
            elif matched_isp_list[i][1][0] != 0 and matched_isp_list[j][1][0] == 0:
                change = 0
                for k in range(len(matched_isp_list[i][1])):
                    if matched_isp_list[i][1][k] != 0 and matched_isp_list[j][1][k] != 0:
                        break
                    if change == 0 and matched_isp_list[i][1][k] != 0 and matched_isp_list[j][1][k] == 0:
                        continue
                    elif change == 0 and matched_isp_list[i][1][k] == 0 and matched_isp_list[j][1][k] == 0:
                        if k == len(matched_isp_list[i][1]) - 1:
                            flag = 1
                        continue
                    elif change == 0 and matched_isp_list[i][1][k] == 0 and matched_isp_list[j][1][k] != 0:
                        change = 1
                        if k == len(matched_isp_list[i][1]) - 1:
                            flag = 1
                    elif change == 1 and matched_isp_list[i][1][k] == 0 and matched_isp_list[j][1][k] == 0:
                        if k == len(matched_isp_list[i][1]) - 1:
                            flag = 1
                        continue
                    elif change == 1 and matched_isp_list[i][1][k] == 0 and matched_isp_list[j][1][k] != 0:
                        if k == len(matched_isp_list[i][1]) - 1:
                            flag = 1
                        continue
                    elif change == 1 and matched_isp_list[i][1][k] != 0 and matched_isp_list[j][1][k] == 0:
                        break

            # 如果可以合并
            if flag == 1:
                # 以i为模板制作一个新的峰簇，用j将其补全
                temp_isp = []
                temp_isp.append(matched_isp_list[i][0].copy())
                temp_isp.append(matched_isp_list[i][1].copy())
                for k in range(len(temp_isp[0])):
                    if matched_isp_list[i][1][k] != 0:
                        continue
                    else:
                        # 同时更新m/z信息和峰强信息
                        temp_isp[0][k] = matched_isp_list[j][0][k]
                        temp_isp[1][k] = matched_isp_list[j][1][k]

                # 计算分数，减去原来的，加上新的
                score1 = matched_score_list[i]
                score2 = matched_score_list[j]
                score3 = score_match_isp(the_isp[index_the[i]], temp_isp, ppm)
                if np.sum(the_lost_list[index_the[i]]) > 0:
                    score3 = p * score3
                global_score = global_score - score1 - score2 + score3
                # 将原来的峰簇标记
                matched_isp_list[i] = 'deleted'
                matched_isp_list[j] = 'deleted'
                matched_score_list[i] = -1
                matched_score_list[j] = -1
                index_the.append(index_the[i])
                index_the[i] = -1
                index_the[j] = -1
                # 记录应删除峰簇的数量
                delete_counter = delete_counter + 2
                # 加入合并后的峰簇
                matched_isp_list.append(temp_isp)
                matched_score_list.append(score3)

            # 无须合并，找匹配分数最大的进行保留
            else:
                score1 = matched_score_list[i]
                score2 = matched_score_list[j]
                if score1 < score2:
                    global_score = global_score - score1
                    matched_isp_list[i] = 'deleted'
                    matched_score_list[i] = -1
                    index_the[i] = -1
                else:
                    global_score = global_score - score2
                    matched_isp_list[j] = 'deleted'
                    matched_score_list[j] = -1
                    index_the[j] = -1
                delete_counter = delete_counter + 1

        else:
            temp_score_list = []
            temp_max_score = 0
            temp_index = -1
            i = 0
            for index_item in item_list:
                score = matched_score_list[index_item]
                temp_score_list.append(score)
                if score > temp_max_score:
                    temp_max_score = score
                    temp_index = i
                i = i + 1

            for j, score in enumerate(temp_score_list):
                if j == temp_index:
                    continue
                global_score = global_score - score
                matched_isp_list[item_list[j]] = 'deleted'
                matched_score_list[item_list[j]] = -1
                index_the[item_list[j]] = -1
                delete_counter = delete_counter + 1

    # 根据标记数量，删除标记峰簇
    for i in range(delete_counter):
        matched_isp_list.remove('deleted')
        index_the.remove(-1)
        matched_score_list.remove(-1)

    # 更新理论谱list
    ret_the_isp_list = []
    ret_the_lost_list = []
    ret_the_z_list = []
    ret_the_HNa_list = []
    for i in index_the:
        ret_the_isp_list.append(the_isp[i])
        ret_the_lost_list.append(the_lost_list[i])
        ret_the_z_list.append(the_z_list[i])
        ret_the_HNa_list.append(the_HNa_list[i])
    score_mass = []  # 标记 得分的单同位素峰
    for i in range(len(matched_isp_list)):
        mass = matched_isp_list[i][0]
        score_mass.append(mass)

    return global_score, matched_isp_list, ret_the_isp_list, matched_score_list, score_mass, ret_the_lost_list, \
           ret_the_z_list, ret_the_HNa_list


def calculate_all_hp_score(the_spectra, dict_list, the_HP, max_int, ppm=10, thresh=0, min_match_num=3, prob=0.5):
    exp_isp = get_exp_isp(dict_list, max_int)
    exp_isp = sort_exp(exp_isp)
    exp_isp_01 = change_sp_format(exp_isp, 0.1)
    all_comp = []
    all_atom = []
    all_comp_mass = []
    all_comp_score = []
    all_comp_match = []
    all_comp_lost = []
    all_comp_z = []
    all_comp_HNa = []
    dict_match_exp = {}
    dict_theo_list = []
    all_score_list = []
    all_score_mass = []
    count = 0
    the_comp_list = the_HP[2]
    the_atom_list = the_HP[3]
    the_mass_list = the_HP[4]
    pass_count = 0
    ttt = time()
    for i in range(len(the_comp_list)):
        count += 1
        print("ID\t" + str(count) + "\t")
        x = the_mass_list[i]
        one_comp = the_comp_list[i]
        one_atom = the_atom_list[i]
        the_isp = the_spectra[0][i]
        lost_comp = the_spectra[1][i]
        Z_list = the_spectra[2][i]
        HNa_list = the_spectra[4][i]
        the_isp_01 = transform_the_01(the_spectra[3][i])
        com_score = compared_score(exp_isp_01, the_isp_01)
        if com_score <= 0:
            # print(com_score)
            pass_count += 1
            continue
        score, match_list, match_theo_list, score_list, score_mass, match_lost_list, match_z_list, match_HNa_list, \
            = match_isotopic(exp_isp, the_isp, ppm, thresh, min_match_num, lost_comp, Z_list,HNa_list, prob)
        if score > 0:
            dict_match_theo = {}
            for i in range(len(match_list)):
                for j in range(len(match_list[i][0])):
                    mass = match_list[i][0][j]
                    exp_int = match_list[i][1][j]
                    the_int = match_theo_list[i][1][j]
                    if mass not in dict_match_exp.keys():
                        dict_match_exp[mass] = exp_int
                    dict_match_theo[mass] = the_int
            dict_theo_list.append(dict_match_theo)
            all_comp.append(one_comp)
            all_atom.append(one_atom)
            all_comp_mass.append(x)
            all_comp_score.append(score)
            all_comp_match.append(match_list)
            all_score_list.append(score_list)
            all_score_mass.append(score_mass)
            all_comp_lost.append(match_lost_list)
            all_comp_z.append(match_z_list)
            all_comp_HNa.append(match_HNa_list)
        #     print(count, one_comp, x, score)
        #     print(len(dict_theo_list))
        # if count % 100 == 0:
        #     print(count)

    print('passed', pass_count)
    # print("ID\t" + str(count + pass_count) + "\t")
    print('ma1', time() - ttt)
    return all_comp, all_atom, all_comp_mass, all_comp_score, dict_match_exp, \
           dict_theo_list, all_score_list, all_score_mass, all_comp_lost, all_comp_z, all_comp_HNa



