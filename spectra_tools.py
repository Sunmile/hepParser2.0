"""
author:houmeijie
author2: wanghui
"""
import colorsys
import os
import random
import subprocess
import platform
import sys
import copy
import numpy as np
import xml.dom.minidom as xml
from brainpy import isotopic_variants
from PyQt5.QtCore import *

from Hp_opt import data_process


dict_atom = {}
dict_atom['H1'] = 1.00783
dict_atom['H2'] = 2.01410
dict_atom['C12'] = 12.00000
dict_atom['C13'] = 13.00335
dict_atom['N14'] = 14.00307
dict_atom['N15'] = 15.00011
dict_atom['O16'] = 15.99492
dict_atom['O17'] = 16.99913
dict_atom['O18'] = 17.99916
dict_atom['S32'] = 31.97207
dict_atom['S33'] = 32.97146
dict_atom['S34'] = 33.96787
dict_atom['S36'] = 35.96708
dict_atom['Na'] = 22.98977

# 用于转换raw文件
class ConvertWorker(QThread):
    sinOut = pyqtSignal(str)

    def __init__(self, parent=None, curPath=None, rawFileName=None):
        super(ConvertWorker, self).__init__(parent)
        self.curPath = curPath
        self.rawFileName = rawFileName

    def run(self):
        os.chdir('msconvert')
        # msconvert ../../../Downloads/rawdata/rawdata -v --zlib --filter "msLevel 1” --filter "peakPicking cwt snr=0.1 peakSpace=0.1 msLevel=1-"
        # command = 'msconvert ' + self.rawFileName + ' -v --zlib --filter "msLevel 1" --filter "peakPicking cwt snr=0.1 peakSpace=0.1 msLevel=1-" -o ' + self.curPath
        command = 'msconvert ' + self.rawFileName + ' -v --zlib --filter "msLevel 1" -o ' + self.curPath
        print(command)
        cmd = subprocess.Popen(command, stdin=subprocess.PIPE, stderr=subprocess.PIPE,
                               stdout=subprocess.PIPE, universal_newlines=True, shell=True, bufsize=1)
        # 实时输出
        while True:
            line = cmd.stdout.readline()
            if line.startswith("converting spectra"):
                proceMess = line.split(" ")[2]
                print(proceMess)
                self.sinOut.emit(proceMess)
            if subprocess.Popen.poll(cmd) == 0:  # 判断子进程是否结束
                break
        os.chdir('../')


# 用于转换mzML文件的TIC图,谱文件借助于pymzml
class mzMLWorker(QThread):
    sinXml = pyqtSignal(list)

    def __init__(self, parent=None, xmlFileName=None):
        super(mzMLWorker, self).__init__(parent)
        self.xmlFileName = xmlFileName

    def run(self):
        xmlInfo = self.loadXML(self.xmlFileName)
        self.sinXml.emit(xmlInfo)

    # 使用minidom解析器打开 XML 文档
    def loadXML(self, fileName):
        DOMTree = xml.parse(fileName)
        root = DOMTree.documentElement

        times, tics, scans = [], [], []
        maxTic, maxScan = 0.0, 0
        num = 0
        scan = 0
        spectrums = root.getElementsByTagName("spectrum")

        for spectrum in spectrums:
            # scan = spectrum.getAttribute("id").split()[2].split('=')[1]
            # msNode, ticNode, timeNode = spectrum.childNodes[3], spectrum.childNodes[13], spectrum.childNodes[15]
            mslevel, tic, startTime = None, None, None
            for i in range(1, len(spectrum.childNodes), 2):
                node = spectrum.childNodes[i]
                if mslevel is None and node.attributes['name'].value == "ms level":
                    mslevel = node.attributes['value'].value
                if tic is None and node.attributes['name'].value == "total ion current":
                    tic = float(node.attributes['value'].value)
                if startTime is None and node.localName == 'scanList':
                    scanNode = node.childNodes[3]
                    for j in range(1, len(scanNode.childNodes), 2):
                        if scanNode.childNodes[j].attributes['name'].value == "scan start time":
                            startTime = float(scanNode.childNodes[j].attributes['value'].value)
                            break

            if mslevel == '1' or mslevel=='2':
                scan += 1
                num += 1
                times.append(startTime)
                tics.append(tic)
                scans.append(scan)
                if tic > maxTic:
                    maxTic = tic
                    maxScan = scan
        print("mzml:", num, maxScan, maxTic, "xmlSpectrums:", len(spectrums))
        return [times, tics, scans, maxScan]


# 用于进行数据处理
class DataWorker(QThread):
    sinID = pyqtSignal(int)
    sinDataInfo = pyqtSignal(list)

    def __init__(self, parent=None, peaks=None, exp_isp=None, max_int=None, total_int=None, the_spectra=None, dict_list=None,
                 the_HP=None, ppm=None,
                 bound_th=None, bound_intensity=None, chb_dA=False, chb_aM=False, chb_aG=False, min_dp=None,max_dp=None):
        super(DataWorker, self).__init__(parent)
        self.peaks = peaks
        self.exp_isp = exp_isp
        self.max_int = max_int
        self.total_int = total_int
        self.the_spectra = the_spectra
        self.dict_list = dict_list
        self.the_HP = the_HP
        self.ppm = ppm
        self.bound_th = bound_th
        self.bound_intensity = bound_intensity
        self.chb_dA = chb_dA
        self.chb_aM = chb_aM
        self.chb_aG = chb_aG
        self.min_dp = min_dp
        self.max_dp = max_dp

    def run(self):
        class redirect:

            def __init__(self, sinID):
                self.sinID = sinID
                self.content = ""

            def write(self, str, content=None):
                if str.startswith("ID"):
                    if int(str.split("\t")[1]) % 100 == 0 or int(str.split("\t")[1]) == 2945:
                        self.sinID.emit(int(str.split("\t")[1]))
                self.content += str
                # else:
                #     self.content += str

        stdout = sys.stdout
        sys.stdout = redirect(self.sinID)

        data_info = data_process(peaks=self.peaks, exp_isp=self.exp_isp, max_int=self.max_int, total_int=self.total_int,
                                 the_spectra=self.the_spectra,
                                 dict_list=self.dict_list, the_HP=self.the_HP,
                                 chb_dA=self.chb_dA, chb_aM=self.chb_aM, chb_aG=self.chb_aG,
                                 min_dp=self.min_dp, max_dp=self.max_dp,
                                 ppm=self.ppm)

        re, sys.stdout = sys.stdout, stdout
        # data_info = [peaks, max_intensity, exp_isp, label_info, candidate_max_num]
        print(re.content)
        self.sinDataInfo.emit(data_info)


# 将谱峰转化为适合画图的数据
def format_peaks(peaks):
    xs, ys = [], []
    for line in peaks:
        xs += line[0], line[0], line[0]
        ys += 0.0, line[1], 0.0
    return xs, ys


# 将谱峰转化为适合画图的数据
def format_peaks_350(peaks):
    xs, ys = [], []
    for line in peaks:
        if line[0] >= 0:
            xs += line[0], line[0], line[0]
            ys += 0.0, line[1], 0.0
    return xs, ys


def format_peaks_2_dim(mzs, intens):
    xs, ys = [], []
    for line in zip(mzs, intens):
        xs += line[0], line[0], line[0]
        ys += 0.0, line[1], 0.0
    return xs, ys


def format_peaks_alone(peaks):
    xs, ys = [], []
    for line in peaks:
        xs.append([line[0], line[0], line[0]])
        ys.append([0.0, line[1], 0.0])
    return xs, ys


# 将同位素簇转为峰列表,[358.123, 567, '-1','3']
# 返回：同为素峰簇的xs,ys,最高同位素簇列表
def get_isotope(peak_dicts, exp_isp):
    isotope_peaks = []
    isotope_max_peaks = []
    for i in range(0, len(exp_isp)):
        mz_list, inten_list = exp_isp[i]
        z = int(1 / (mz_list[1] - mz_list[0]) + 0.5)  # 四舍五入
        max_index, max_inten = 0, 0
        for j in range(0, 5):
            if inten_list[j] != 0:
                isotope_peaks.append([mz_list[j], peak_dicts[mz_list[j]]])
                if max_inten < inten_list[j]:
                    max_index = j
                    max_inten = inten_list[j]
        isotope_max_peaks.append([mz_list[max_index], peak_dicts[mz_list[max_index]], z, j + 1])

    isotope_peaks.sort(key=lambda x: x[0])

    isotope_xs, isotope_ys = format_peaks(isotope_peaks)
    return isotope_xs, isotope_ys, isotope_max_peaks


# 转换标注的信息
# label是约定好的六元组，[[mass, Z, HNa,component, lose , score],[mass, Z HNa, component, lose , score]]
# dict_mass_comp,key:理论中性质量，value: 理论结构分子组成
# dict_mass_flag,key:理论中性质量，value: flag 数组， flag[0]: 实验谱匹配m/z,OR 理论中性mass-1, flag[1]: 电荷数，OR 0
# dict_mass_family,key:理论中性质量，value: 衍生 label 数组 index

def get_all_struct(label_info):
    right_comp = []
    label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order, df_conf = label_info
    ascend_key = sorted(key_with_order)
    for mass in ascend_key:
        tmp_index = dict_mass_family[mass]
        comp = dict_mass_comp[mass]

        per_struct = ""
        # for ele in comp:
        #     per_struct += str(ele) + ","
        # per_struct = "[" + per_struct[:-1] + "]"
        if comp[0] == 1:
            per_struct += 'dA'
        per_struct += '['+ str(comp[1])+','+str(comp[2])+','+str(comp[3])+','+str(comp[4]) + ']'
        if comp[5] == 1:
            per_struct += 'aG'
        if comp[6] == 1:
            per_struct += 'aM'
        right_comp.append([per_struct, '{:.2f}'.format(label[tmp_index[0]][5])])

    return right_comp, key_with_order

def get_molecular_mass(atom_list):
    Mass = 0
    atom_mass = [dict_atom['H1'], dict_atom['C12'], dict_atom['N14'],
                 dict_atom['O16'], dict_atom['S32'], dict_atom['Na']]
    for i in range(len(atom_list)):
        Mass += atom_mass[i] * atom_list[i]
    return np.round(Mass, 5)

def get_theory_mz(z, hna, comp, lost,is_HNa='H'):
    comp_atoms = [
        [8, 6, 0, 6, 0],  # 不饱和糖醛酸 HexA
        [10, 6, 0, 7, 0],  # 饱和糖醛酸  GlcA
        [13, 6, 1, 5, 0],  # 葡萄糖胺   GlcN
        [4, 2, 0, 2, 0],  # 乙酰基    Ac
        [2, 0, 0, 4, 1],  # 硫酸基    SO3
        [11, 6, 1, 4, 0],  # 内醚糖    Levoglucosan
        [12, 6, 0, 5, 0]  # 甘露糖(dicp jin)    Man
    ]
    comp_atoms = np.array(comp_atoms)
    tmp_atoms = np.array([2, 0, 0, 1, 0])  # 带一个水
    H2O = np.array([2, 0, 0, 1, 0])  # 任两个基团结合，都脱一个水
    for i in range(len(comp)):
        tmp_atoms = tmp_atoms + comp_atoms[i] * comp[i] - H2O * comp[i]
    lost_list = [[0, 0, 0, 3, 1],  # 丢失SO3
                 [1, 0, 1, 0, 0],  # 丢失NH
                 [1, 0, 1, 3, 1],  # 丢失NHSO3
                 [0, 1, 0, 2, 0]]  # 丢失COO
    lost_list = np.array(lost_list)
    for i in range(len(lost)):
        tmp_atoms = tmp_atoms - lost_list[i] * lost[i]
    tmp_com = tmp_atoms
    mass = get_molecular_mass(tmp_atoms)
    tmp_lwh = {}
    tmp_lwh['H'] = tmp_com[0]
    tmp_lwh['C'] = tmp_com[1]
    tmp_lwh['N'] = tmp_com[2]
    tmp_lwh['O'] = tmp_com[3]
    tmp_lwh['S'] = tmp_com[4]
    tmp_distribution = isotopic_variants(tmp_lwh, npeaks=5, charge=-z)
    tmp_peak = tmp_distribution[0]
    add_mass=0
    if is_HNa =='HNa':
        add_mass = dict_atom['Na'] - 1
    elif is_HNa =='HNH4':
        add_mass = dict_atom['N14'] + 4 * dict_atom['H1'] - 1
    diff_mass = hna[1] * (add_mass)  # 计算与只脱H的mass的差距
    diff_mass = diff_mass / np.abs(z)
    mz = np.round(tmp_peak.mz + diff_mass, 4)
    return mass, mz


def get_labels(label_info, peak_dict, lose_ion, new_key_with_order, is_HNa='H'):
    xy = []
    mass_label = {}  # (mz,inten):[mz,Z,compoment,lose,score],1和0分别是匹配上了原始结构，没匹配上原始结构
    mass_struct_tips = {}
    struct_id = {}
    comp_copy = []
    key_id = {}
    right_comp = []
    new_labels = []
    label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order, df_conf = label_info
    num = 0

    for mass in new_key_with_order:
        tmp_index = dict_mass_family[mass]
        comp_copy.append([mass, label[tmp_index[0]][5]])

    comp_copy.sort(key=lambda x: x[1], reverse=True)
    for i, element in enumerate(comp_copy):
        key_id[element[0]] = i + 1

    for mass in new_key_with_order:
        num += 1
        tmp_index = dict_mass_family[mass]
        comp = dict_mass_comp[mass]
        p_adj = df_conf[df_conf.comp==str(comp)]['log_p'].values[0]
        per_struct = ""
        # for ele in comp:
        #     per_struct += str(ele) + ","
        if comp[0] == 1:
            per_struct += 'dA'
        per_struct += '['+ str(comp[1])+','+str(comp[2])+','+str(comp[3])+','+str(comp[4]) + ']'
        if comp[5] == 1:
            per_struct += 'aG'
        if comp[6] == 1:
            per_struct += 'aM'
        # per_struct = "[" + per_struct[:-1] + "]"
        right_comp.append([per_struct, '{:.2f}'.format(label[tmp_index[0]][5]),np.round(p_adj,2), num, key_id[mass]])
        struct_id[per_struct] = num
        for id in tmp_index:
            mz, z, hna, _, lose, score, th = label[id]  # _就是comp
            xy.append((mz, peak_dict[mz]))
            if (mz, peak_dict[mz]) not in mass_label.keys():
                mass_label[(mz, peak_dict[mz])] = []  # 分别是分子式序号和颜色序号
            hna = np.array(hna)
            mass_label[(mz, peak_dict[mz])].append([comp, lose, z, hna, score, th, num, key_id[mass]])

    # 标注的峰
    xy = list(set(xy))
    xy.sort(key=lambda x: x[0])

    # 标注的tips
    for key in mass_label.keys():
        struct_info, structs, scores = format_struct_2(mass_label[key], lose_ion, is_HNa)
        mass_struct_tips[key] = "$m/z : " + str(key[0]) + "$ \n$intensity : " + str(key[1]) + "$\n" + struct_info

        # 转换label为5个list
    mass_list, int_list, z_list, hna_list, the_mz_list, comp_list, the_mass_list, lose_list, score_list, th_list = \
        [], [], [], [], [], [], [], [], [], []

    label_copy = copy.deepcopy(label)
    label_copy.sort(key=lambda x: x[0])
    for line in label_copy:
        per_struct = ""
        if line[3][0] == 1:
            per_struct += 'dA'
        per_struct += '['+ str(line[3][1])+','+str(line[3][2])+','+str(line[3][3])+','+str(line[3][4]) + ']'
        if line[3][5] == 1:
            per_struct += 'aG'
        if line[3][6] == 1:
            per_struct += 'aM'
        new_labels.append([line[0], line[1], line[2][0],line[2][1], line[6], per_struct, line[4]])
        mass_list.append(line[0])
        int_list.append(peak_dict[line[0]])
        z_list.append(line[1])
        hna_list.append(line[2])
        tmp_mass,tmp_mz = get_theory_mz(line[1],line[2],line[3],line[4], is_HNa=is_HNa)
        the_mz_list.append(tmp_mz)
        the_mass_list.append(tmp_mass)
        comp_list.append(per_struct)
        lose_list.append(line[4])
        # score_list.append(line[4])
        th_list.append(line[6])
    hna_list = np.array(hna_list)
    return xy, mass_label, mass_struct_tips, right_comp, struct_id, \
           (mass_list,int_list, z_list,hna_list,the_mz_list, th_list, comp_list, lose_list, the_mass_list), new_labels


# 寻找衍生峰的mz、inten
def get_family(label_info, peak_dict):
    label, dict_mass_comp, dict_mass_flag, dict_mass_family, key_with_order, df_conf = label_info
    mass_family = {}
    for key in dict_mass_family.keys():
        mz = dict_mass_flag[key][0]
        index_list = dict_mass_family[key]
        mz_list, inten_list, comp_list, lose_list, z_list,hna_list, score_list = [], [], [], [], [], [],[]
        for index in index_list:
            f_mz, Z, hna,component, lose, score, th = label[index]
            mz_list.append(f_mz)
            inten_list.append(peak_dict[f_mz])
            comp_list.append(component)
            lose_list.append(lose)
            z_list.append(Z)
            hna_list.append(hna)
            score_list.append(score)
        mass_family[key] = [mz_list, inten_list, comp_list, lose_list, z_list, hna_list, score_list]
    return mass_family


# 寻找最近的数据点
def find_closed_key(comps_xys, x, xys):
    min_dis = 1000
    min_index = 0
    for i in range(len(comps_xys)):
        if abs(comps_xys[i][0] - x) < min_dis:
            min_dis = abs(comps_xys[i][0] - x)
            min_index = i
    return xys[min_index]


# 生成分子式的坐标
def generate_annotation_xy(min_x, max_x, n, y):
    annotation_xys = []
    y = 1.1 * y
    delta = (max_x - min_x) / (n - 1)
    for i in range(0, n):
        annotation_xys.append((min_x + i * delta, y))
    return annotation_xys


# 利用上下空间区分开分子式的标注
def generate_up_down_xy(xys, up_x_percentage, up_y_percentage):
    standard_xy = []
    if len(xys) == 0:
        return standard_xy
    if len(xys) == 1:
        standard_xy.append([xys[0]])
        return standard_xy

    diff_x_list, diff_y_list = [], []
    xys.sort(key=lambda x: x[0])
    for i in range(1, len(xys)):
        diff_x_list.append(xys[i][0] - xys[i - 1][0])
    diff_x_list.sort()
    if up_x_percentage >= 1:
        x_delta = diff_x_list[-1]
    else:
        # print(diff_x_list)
        x_delta = diff_x_list[int(len(diff_x_list) * up_x_percentage)]

    tmp_list = [xys[0]]
    for i in range(1, len(xys)):
        if xys[i][0] - xys[i - 1][0] > x_delta or (xys[i][0] - xys[i - 1][0] <= x_delta and (
                xys[i][1] / xys[i - 1][1] > (1 + up_y_percentage)
                or xys[i - 1][1] / xys[i][1] > (1 + up_y_percentage))):
            standard_xy.append(tmp_list)
            tmp_list = [xys[i]]
        else:
            tmp_list.append(xys[i])
    standard_xy.append(tmp_list)
    return standard_xy


def format_struct(label_list, lose_ion):
    struct_info = ""
    structs = []
    scores = []
    for label in label_list:
        per_struct = ""
        for ele in label[0]:
            per_struct += str(ele) + ","
        per_struct = "[" + per_struct[:-1] + "]"
        per_lost = ""
        for i in range(0, 4):
            if label[1][i] == 0:
                continue
            elif label[1][i] == 1:
                per_lost += "-" + lose_ion[i]
            else:
                per_lost += "-" + str(label[1][i]) + lose_ion[i]
        per_struct += per_lost + "<sup>" + str(label[2]) + "-</sup>"
        structs.append(per_struct)
        struct_info += per_struct
        scores.append(round(label[-1], 5))
    return struct_info, structs, scores


def format_struct_2(label_list, lose_ion,is_HNa):
    struct_info = ""
    structs = []
    scores = []
    for i in range(len(label_list) - 1, -1, -1):
        label = label_list[i]
        # for label in label_list:
        per_struct = ""
        # for ele in label[0]:
        #     per_struct += str(ele) + ","
        # per_struct = "[" + per_struct[:-1] + "]"
        if label[0][0] == 1:
            per_struct += 'dA'
        per_struct += '['+ str(label[0][1])+','+str(label[0][2])+','+str(label[0][3])+','+str(label[0][4]) + ']'
        if label[0][5] == 1:
            per_struct += '^{aG}'
        if label[0][6] == 1:
            per_struct += '_{aM}'
        per_lost = ""
        for i in range(0, 4):
            # print(label[1])
            if label[1][i] == 0:
                continue
            elif label[1][i] == 1:
                per_lost += "-" + lose_ion[i]
            else:
                per_lost += "-" + str(label[1][i]) + lose_ion[i]
        if per_lost=='' and label[0][5]==1:
            per_struct = per_struct[:-1] +str(label[2])+"-}"
        else:
            per_struct += per_lost + "^{" + str(label[2]) + "-}"
        tmp_hna = label[3]
        per_struct += '-'+str(tmp_hna[0])+'H+'+str(tmp_hna[1])+is_HNa[1:]
        structs.append(per_struct)
        struct_info += '$' + str(label[-2]) + " : " + per_struct + '$ \n'
        scores.append(round(label[4], 5))
    struct_info = struct_info + " "
    return struct_info[:-2], structs, scores


def get_first_max_num(score_list, key_with_order):
    select_id = []
    if len(score_list)<1:
        return select_id, 0
    max_score_id = np.argmax(score_list)
    max_score = np.max(score_list)

    for i in range(len(score_list)):
        if i <= max_score_id:
            id = key_with_order[i] // 10000
            select_id.append(id)
    print("better:", max_score)
    return select_id, max_score


def format_loss(lose_list, lose_ion):
    per_lost = ""
    for j in range(0, 4):
        if lose_list[j] == 0:
            continue
        elif lose_list[j] == 1:
            per_lost += "-" + lose_ion[j]
        else:
            per_lost += "-" + str(lose_list[j]) + lose_ion[j]
    return "$" + per_lost + "$"


def format_comp(comp_list, z):
    comp_info = ""
    if comp_list[0] == 1:
        comp_info += 'dA'
    comp_info += '[' + str(comp_list[1]) + ',' + str(comp_list[2]) + ',' + str(comp_list[3]) + ',' + str(
        comp_list[4]) + ']'
    if comp_list[5] == 1:
        comp_info += '^{aG'+str(z)+'-}'
    if comp_list[6] == 1:
        comp_info += '_{aM}'
    # for i in range(len(comp_list)):
    #     comp_info += str(comp_list[i]) + ","
    if comp_list[5] == 1:
        comp_info = "$" + comp_info + "$"
    else:
        comp_info = "$" + comp_info + "^{" + str(z) + "-}$"
    return comp_info


def get_n_hls_colors(num):
    hls_colors = []
    i = 0
    step = 360.0 / num
    while i < 360:
        h = i
        s = 90 + random.random() * 10
        l = 50 + random.random() * 10
        _hlsc = [h / 360.0, l / 100.0, s / 100.0]
        hls_colors.append(_hlsc)
        i += step

    return hls_colors


def format_color(color):
    digit = list(map(str, range(10))) + list("ABCDEF")
    string = '#'
    for i in range(0, 3):
        a1 = int(color[i] * 255) // 16
        a2 = int(color[i] * 255) % 16
        string += digit[a1] + digit[a2]
    return string


def dataClosed():
    if platform.system() == "Windows":
        os.system("rmdir /s/q mzML")
        os.system("mkdir mzML")
    else:
        pass

    # 颜色1：黑色，原始峰   struct_info, structs, scores = format_struct(self.massToLabels[key], self.lose_ion)
    # 颜色2：绿色, 同位素簇 label.sticks
    # 颜色3：绿色，最高的同位素标记电荷 label.sticks
    # 颜色4: 蓝色, 最高的同位素点击时提示 分子式，但并不直接写上   hoverinfo
