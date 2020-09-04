from PyQt5.QtCore import *
import xml.dom.minidom as xml
import binascii
import zlib
import time
import numpy as np
import pickle as pk


class Peak:
    def __init__(self, mzValue, intensValue):
        self.mzValue = mzValue
        self.intensValue = intensValue


class Spectrum:
    def __init__(self, scan, mslevel, pepmass):
        self.scan = scan
        self.mslevel = mslevel
        self.pepmass = pepmass


# 用于转换mzML文件的TIC图,谱文件借助于pymzml
class mzMLWorker(QThread):
    sinXml = pyqtSignal(list)

    def __init__(self, parent=None, xmlFileName=None):
        super(mzMLWorker, self).__init__(parent)
        self.xmlFileName = xmlFileName
        self.rootPath = "./"

    def run(self):
        xmlInfo = self.loadXML()
        self.sinXml.emit(xmlInfo)

    def loadXML(self):
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

        print("mzml:", num)
        t1 = time.time()
        print("mzML to Mgf:", t1 - t)
        draw_data = self.get_dict(peaks, tics)
        print("draw iter", time.time() - t1)

        return [times, tics, scans, maxScan, draw_data]

    # 先对数据进行zlib压缩，再进行Base64加密，解密时先解Base64，再解zlib
    def decodeBase64AndDecompressZlib(self, data):
        time0 = time.time()
        dec_data = binascii.a2b_base64(data)
        dec_data = zlib.decompress(dec_data)
        ret_data = np.frombuffer(dec_data, np.float64)

        return ret_data

    def save_mgf(self, scan, mz_list, inten_list):
        print(scan)
        with open(self.rootPath + "/" + str(scan), 'w', newline='') as f:
            for i in range(0, len(mz_list)):
                f.write(str(mz_list[i]) + " " + str(inten_list[i]) + "\n")

    def get_dict(self, peaks, tic):
        THRESHOLD = 0.001
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
        for i in range(maxmass + 1):
            z.append([])
            x.append(i)

        for i, peak_list in enumerate(peaks):
            max_peak = 0
            for tpeak in peak_list[1]:
                if tpeak > max_peak:
                    max_peak = tpeak

            max_peak = max_peak * THRESHOLD
            y.append(i)
            for k in range(len(z)):
                z[k].append(0)
            for j, tpeak in enumerate(peak_list[1]):
                if tpeak < max_peak:
                    continue

                t_mass0 = int(round(peak_list[0][j]))
                z[t_mass0][-1] += int(tpeak)

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

        ret_z = self.smooth(ret_z)
        with open('data/dict.pk', 'wb') as f:
            pk.dump([ret_x, y, ret_z, tic], f)
        return [ret_x, y, ret_z, tic]

    def smooth(self, z_data):
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
    ret_data = np.frombuffer(dec_data, np.float64)

    return ret_data


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
        for (mz, inten) in zip(mz_list, inten_list):
            new_mz = round(mz, 4)
            if new_mz not in mz_inten.keys():
                mz_inten[new_mz] = 0
            mz_inten[new_mz] += inten
    for key in mz_inten.keys():
        peaks.append([key, mz_inten[key]])
        max_inten = max(max_inten, mz_inten[key])
    peaks.sort(key=lambda x: x[0])
    t1 = time.time()
    print("get_merged_peaks", t1 - t)
    return peaks, max_inten

# parser = XmlParser("")
# infos = parser.loadXML()