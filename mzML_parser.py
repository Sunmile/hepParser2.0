from PyQt5.QtCore import *
import xml.dom.minidom as xml
import binascii
import zlib
import time
import numpy as np
import pickle as pk
import threading


class Peak:
    def __init__(self, mzValue, intensValue):
        self.mzValue = mzValue
        self.intensValue = intensValue


class Spectrum:
    def __init__(self, scan, mslevel, pepmass):
        self.scan = scan
        self.mslevel = mslevel
        self.pepmass = pepmass


def add_num_thread(z, peak_list, begin, done):
    # threading.Lock.acquire()
    THRESHOLD = 0.01
    for i in range(len(peak_list)):
        max_peak = 0
        if len(peak_list[i][1]) == 0:
            continue
        max_peak = np.max(np.array(peak_list[i][1]))
        max_peak = max_peak * THRESHOLD
        for j, tpeak in enumerate(peak_list[i][1]):
            if tpeak < max_peak:
                continue
            t_mass0 = int(round(peak_list[i][0][j]))
            z[t_mass0][i + begin] += int(tpeak)
    done[0] += len(peak_list)


# 用于转换mzML文件的TIC图,谱文件借助于pymzml
class mzMLWorker(QThread):
    sinXml = pyqtSignal(list)
    step = pyqtSignal(int)

    def __init__(self, parent=None, xmlFileName=None):
        super(mzMLWorker, self).__init__(parent)
        self.xmlFileName = xmlFileName
        self.rootPath = "./"

    def run(self):
        xmlInfo = self.loadXML()
        self.sinXml.emit(xmlInfo)

    def loadXML(self):
        self.step.emit(0)
        DOMTree = xml.parse(self.xmlFileName)
        root = DOMTree.documentElement

        times, tics, scans = [], [], []
        maxTic, maxScan = 0.0, 0
        num = 0
        spectrums = root.getElementsByTagName("spectrum")
        peaks = []
        t = time.time()
        for spectrum in spectrums:
            scan = spectrum.getAttribute("id").split()[2].split('=')[1]
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

            binary_nodes = spectrum.getElementsByTagName("binary")

            if mslevel == '1':
                num += 1
                times.append(startTime)
                tics.append(tic)
                scans.append("scan=" + scan)
                if tic > maxTic:
                    maxTic = tic
                    maxScan = int(scan)

                if len(binary_nodes) > 0 and binary_nodes[0].firstChild is not None \
                        and binary_nodes[1].firstChild is not None:
                    mz_data = binary_nodes[0].firstChild.data
                    i_data = binary_nodes[1].firstChild.data
                    mz_list = self.decodeBase64AndDecompressZlib(mz_data)
                    inten_list = self.decodeBase64AndDecompressZlib(i_data)
                    peaks.append([mz_list, inten_list])
                else:
                    peaks.append([[], []])
                    mz_list = []
                    inten_list = []
                # self.save_mgf(scan, mz_list, inten_list)
        self.step.emit(1)
        print("mzml:", num)
        t1 = time.time()
        print("mzML to Mgf:", t1 - t)
        draw_data = self.get_dict(peaks, tics)
        self.step.emit(3)
        print("draw iter", time.time() - t1)
        with open('data/peaks.pk', 'wb') as f:
            pk.dump(peaks, f)

        return [times, tics, scans, maxScan, draw_data]

    # 先对数据进行zlib压缩，再进行Base64加密，解密时先解Base64，再解zlib
    def decodeBase64AndDecompressZlib(self, data):
        time0 = time.time()
        dec_data = binascii.a2b_base64(data)
        dec_data = zlib.decompress(dec_data)
        ret_data = np.frombuffer(dec_data, np.float64)# np.float64

        return ret_data

    def save_mgf(self, scan, mz_list, inten_list):
        print(scan)
        with open(self.rootPath + "/" + str(scan), 'w', newline='') as f:
            for i in range(0, len(mz_list)):
                f.write(str(mz_list[i]) + " " + str(inten_list[i]) + "\n")

    def get_dict(self, peaks, tic):
        THRESHOLDMAXPEAK = 0.05
        maxmass = 0
        for i, peak_list in enumerate(peaks):
            if len(peak_list) < 1:
                continue
            if len(peak_list[0]) < 1:
                continue
            if peak_list[0][-1] > maxmass:
                maxmass = peak_list[0][-1]
        x = []
        y = []
        z = []
        maxmass = int(round(maxmass))
        L = len(peaks)
        for i in range(maxmass + 1):
            z.append([])
            for j in range(L):
                z[-1].append(0)
            x.append(i)

        for i in range(L):
            y.append(i)

        thread_list = []
        batch = 100
        done = [0]
        for i in range(int(len(peaks) / batch) + 1):
            end = min(len(peaks), (i + 1) * batch)
            tthread = threading.Thread(target=add_num_thread, args=(z, peaks[i * batch:end], i * batch, done,))
            thread_list.append(tthread)

        for thread in thread_list:
            thread.start()

        while (done[0] != len(peaks)):
            time.sleep(0.1)

        ret_z = []
        ret_x = []
        max_list = []
        all_max = 0
        for mass_list in z:
            tmax = 0
            for peak in mass_list:
                if peak > tmax:
                    tmax = peak
            if tmax > all_max:
                all_max = tmax
            max_list.append(tmax)

        all_max = all_max * THRESHOLDMAXPEAK
        for i, tmax in enumerate(max_list):
            if tmax < all_max:
                continue
            ret_z.append(z[i])
            ret_x.append(x[i])

        ret_z = smooth(ret_z)
        self.step.emit(2)
        # with open('data/dict.pk', 'wb') as f:
        #    pk.dump([ret_x, y, ret_z, tic], f)
        return [ret_x, y, ret_z, tic]


def smooth(z_data):
    ret = []
    smooth_len = 10
    for z_list in z_data:
        t_list = []
        for i in range(len(z_list)):
            left = max(0, i - int(smooth_len / 2))
            right = min(len(z_list), i + int(smooth_len / 2))

            def ave(data_list):
                L = len(data_list)
                all_d = 0
                for item in data_list:
                    all_d += item
                return int(all_d / L)

            t_list.append(ave(z_list[left:right]))
        ret.append(t_list)

    return ret


def decodeBase64AndDecompressZlib(data):
    dec_data = binascii.a2b_base64(data)
    dec_data = zlib.decompress(dec_data)
    ret_data = np.frombuffer(dec_data, np.float64) # np.float64

    return ret_data


def save_mgf(root_path, scan, mz_list, inten_list):
    print(scan)
    with open(root_path + str(scan), 'w', newline='') as f:
        for i in range(0, len(mz_list)):
            f.write(str(mz_list[i]) + " " + str(inten_list[i]) + "\n")


def get_peaks_by_id(xmlPath, spectrumID):
    DOMTree = xml.parse(xmlPath)
    root = DOMTree.documentElement
    spectrums = root.getElementsByTagName("spectrum")
    spectrum = spectrums[spectrumID - 1]
    binary_nodes = spectrum.getElementsByTagName("binary")
    peaks = []
    maxIntensiy = 0.0
    if len(binary_nodes) > 0 and binary_nodes[0].firstChild is not None \
            and binary_nodes[1].firstChild is not None:
        mz_data = binary_nodes[0].firstChild.data
        i_data = binary_nodes[1].firstChild.data
        mz_list = decodeBase64AndDecompressZlib(mz_data)
        inten_list = decodeBase64AndDecompressZlib(i_data)
    else:
        mz_list = []
        inten_list = []
    for (mz, inten) in zip(mz_list, inten_list):
        peaks.append([mz, inten])
        maxIntensiy = max(maxIntensiy, inten)
    return peaks, maxIntensiy


def get_merged_peaks(xml_path, start_scan, end_scan):
    DOMTree = xml.parse(xml_path)
    root = DOMTree.documentElement
    spectrums = root.getElementsByTagName("spectrum")
    mz_inten = {}
    peaks = []
    max_inten = 0
    t = time.time()
    for i in range(start_scan, end_scan + 1):
        t0 = time.time()
        spectrum = spectrums[i]
        binary_nodes = spectrum.getElementsByTagName("binary")
        if len(binary_nodes) > 0 and binary_nodes[0].firstChild is not None \
                and binary_nodes[1].firstChild is not None:
            mz_data = binary_nodes[0].firstChild.data
            i_data = binary_nodes[1].firstChild.data
            mz_list = decodeBase64AndDecompressZlib(mz_data)
            inten_list = decodeBase64AndDecompressZlib(i_data)
        else:
            mz_list = []
            inten_list = []
        t1 = time.time()
        print("parser_range_peaks", t1 - t0)
        mz_list = np.round(np.array(mz_list), 4)
        for (new_mz, inten) in zip(mz_list, inten_list):
            if new_mz not in mz_inten.keys():
                mz_inten[new_mz] = 0
            mz_inten[new_mz] += inten
        t2 = time.time()
        print("merge", t2 - t1)

    for key in mz_inten.keys():
        peaks.append([key, mz_inten[key]])
        max_inten = max(max_inten, mz_inten[key])
    peaks.sort(key=lambda x: x[0])
    t1 = time.time()
    print("get_merged_peaks", t1 - t)
    return peaks, max_inten


def get_mass_data(file_name, aim_list, digits):
    with open('data/peaks.pk', 'rb') as f:
        peaks = pk.load(f)
    ret_x = list(np.round(np.array(aim_list) * 1.0, digits))
    L = len(peaks)
    ret_y = []
    ret_z = []
    for i in range(L):
        ret_y.append(i)

    dis = 10 ** (-digits)
    for mass in ret_x:
        value = []
        for peak_list in peaks:
            if len(peak_list) == 0:
                value.append(0)
                continue

            if len(peak_list[0]) == 0:
                value.append(0)
                continue

            tmp = 0
            aim_begin = mass - dis * 2
            aim_end = mass + dis * 2
            begin_i = get_index_from_mass_list(peak_list[0], aim_begin)
            end_i = get_index_from_mass_list(peak_list[0], aim_end)
            assert (end_i >= begin_i)
            i = begin_i
            while (i <= end_i):
                a = round(peak_list[0][i], digits)
                b = round(mass, digits)
                if a == b:
                    tmp += int(peak_list[1][i])

                elif a > b:
                    break
                i = i + 1
            value.append(tmp)

        ret_z.append(value)

    ret_z = smooth(ret_z)
    return ret_x, ret_y, ret_z


def get_index_from_mass_list(mass_list, aimnum):
    left = 0
    right = len(mass_list) - 1
    while ((right - left) > 1):
        mid = int((left + right) / 2)
        if mass_list[mid] > aimnum:
            right = mid
        else:
            left = mid

    return left


def loadXML(xmlFileName):
    DOMTree = xml.parse(xmlFileName)
    root = DOMTree.documentElement

    times, tics, scans = [], [], []
    maxTic, maxScan = 0.0, 0
    num = 0
    spectrums = root.getElementsByTagName("spectrum")
    peaks = []
    t = time.time()
    for spectrum in spectrums:
        scan = spectrum.getAttribute("id").split()[2].split('=')[1]
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

        binary_nodes = spectrum.getElementsByTagName("binary")

        if mslevel == '1':
            num += 1
            times.append(startTime)
            tics.append(tic)
            scans.append("scan=" + scan)
            if tic > maxTic:
                maxTic = tic
                maxScan = int(scan)

            if len(binary_nodes) > 0 and binary_nodes[0].firstChild is not None \
                    and binary_nodes[1].firstChild is not None:
                mz_data = binary_nodes[0].firstChild.data
                i_data = binary_nodes[1].firstChild.data
                mz_list = decodeBase64AndDecompressZlib(mz_data)
                inten_list = decodeBase64AndDecompressZlib(i_data)
                peaks.append([mz_list, inten_list])
            else:
                peaks.append([[], []])
                mz_list = []
                inten_list = []
            save_mgf("hep/", scan, mz_list, inten_list)
    print("mzml:", num)
    t1 = time.time()
    print("mzML to Mgf:", t1 - t)


# loadXML("/Users/hou/Downloads/20201117_HEP_F2-F2_CHAI_3uL_NEG IDA_001.mzML")
# parser = XmlParser("")
# infos = parser.loadXML()
