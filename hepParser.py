"""
author:houmeijie
"""
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

from matplotlib.backend_bases import *
from matplotlib.backends.backend_qt5agg import (FigureCanvas, NavigationToolbar2QT as NavigationToolbar)

from candidate_dialog import QmyDialogSize, QAnalyseBar, pbar_str, DPstepBar
from spectra_tools import *
from Hp_opt import get_fit_pk, get_comp_pk, save_file, get_filter_MZ
from PIL import Image
from mzML_parser import *
from draw_3D import *
from qss import *

import matplotlib as mpl
import numpy as np
import pandas as pd
import platform
import mplcursors
import time
import xlwt

if platform.system() != "Windows":
    # 分析界面左侧按钮
    label_btn_str = """QPushButton {
                        background-color: #7b94b1;
                        color :white;
                        padding: 1px 1px;
                        border-radius:10px;
                        margin: 0px;
                    }
                    """
    # 分析界面左侧按钮被按下时
    push_label_btn_str = """QPushButton {
                            background-color: #7b94b1;
                            color :blue;
                            padding: 1px 1px;
                            border-radius:10px;
                            margin: 0px;
                        }
                        """
    # 原始谱界面的 scan编号样式 被按下的样式
    push_scan_label_str = """
                            QPushButton{ 
                                border-radius: 16px;
                                background-color:#3265f6;
                                border:2px solid #3265f6;
                                color: white;
                            }
                            """
    # tic图界面左侧的scan的样式
    tmp_scan_label_str = """
                        QPushButton{
                            border-radius: 16px;
                            padding: 5px 12px;
                            background-color: white;
                            color: #3265f6;
                            border:2px solid #3265f6;
                        }
                        QPushButton:hover {
                            background-color: #3265f6;
                            color:white;
                        }
                        QPushButton:pressed {
                            background-color: #3265f6;
                            color:white;
                        }
                        """
else:
    label_btn_str = """QPushButton {
                            background-color: #7b94b1;
                            color :white;
                            padding: 1px 1px;
                            border-radius:9px;
                            margin: 0px;
                        }
                        """
    # 分析界面左侧按钮被按下时
    push_label_btn_str = """QPushButton {
                                background-color: #7b94b1;
                                color :blue;
                                padding: 1px 1px;
                                border-radius:9px;
                                margin: 0px;
                            }
                            """
    # 原始谱界面的 scan编号样式 被按下的样式
    push_scan_label_str = """
                            QPushButton{ 
                                border-radius: 15px;
                                background-color:#3265f6;
                                border:2px solid #3265f6;
                                color: white;
                            }
                            """
    tmp_scan_label_str = """
                            QPushButton{
                                border-radius: 15px;
                                padding: 5px 12px;
                                background-color: white;
                                color: #3265f6;
                                border:2px solid #3265f6;
                            }
                            QPushButton:hover {
                                background-color: #3265f6;
                                color:white;
                            }
                            QPushButton:pressed {
                                background-color: #3265f6;
                                color:white;
                            }
                            """


class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self._init_config()

        self._init_data()

        self._init_ui()

    def _init_config(self):
        self.init_finish = False
        self.figFlag = [0]  # tic：1  labels:2  famliy:3
        self.total_comp = 5157
        self.ppm = 20
        self.bound_th = 0.001
        self.bound_intensity = 300
        self.x_space = 12
        self.y_space = 20
        self.anno_text_size = 14
        self.opacity = [False, False, False, False]
        self.scan_com_button = {}
        self.chk_list = []
        self.lose_ion = ["HSO_3", 'NH_2', 'NH_{2}SO_3', 'COOH']
        if platform.system() != "Windows":
            self.x_space = 15
            self.y_space = 36
            self.anno_text_size = 13
            print("OS", platform.system())

    def _init_data(self):
        s = time.time()
        self.DATA0 = 'data/plot0_005_tic.pk'

        self.data_3d_file_name = 'data/plot0_005_smooth_tic.pk'
        self.fit_list = get_fit_pk('data/Isotope_dist_fit.pk')
        self.the_HP = get_comp_pk('data/HP_20.pk')  # 枚举的理论肝素结构
        self.peak_dict = {}  # 保存实验谱mz->absolute intensity

        # 用于保存excel
        self.wb = xlwt.Workbook()
        self.ws = self.wb.add_sheet('test')

        e = time.time()
        print("load pk:", e - s)

    def _init_ui(self):

        # 获取总的屏幕大小
        screen = QDesktopWidget().screenGeometry()
        self.setWindowTitle("hepParser")
        self.window().setObjectName("window")
        self.resize(screen.width() * 0.9, screen.height() * 0.9)

        # 设置窗体在屏幕中间
        windowSize = self.geometry()
        self.move((screen.width() - windowSize.width()) / 2, (screen.height() - windowSize.height()) / 8)
        self.setWindowIcon(QIcon('icon/analysing.svg'))

        # 设置外观属性
        mpl.rcParams['font.sans-serif'] = ['SimHei', 'STSong', 'STHeiti', 'Songti SC', 'PingFang HK']  # 汉字字体
        mpl.rcParams['font.size'] = 14  # 字体大小
        mpl.rcParams['axes.unicode_minus'] = False  # 正常显示符号

        # 菜单栏
        self._initMenu()

        # 工具栏
        self._initToolBar()

        # 窗体主要内容
        self._initcentralWidget()

        # 右键菜单栏
        self._init_right_pick_menu()

    def _initcentralWidget(self):

        # 左侧第1页
        self.centralWidget = QWidget(self)

        """左侧边栏的子面板1"""
        # 左侧第1页式tic页面，已被删除

        # 左侧第2页
        self.right_label_page = QWidget()
        self.right_label_page.setGeometry(QRect(0, 0, 356, 401))
        self.right_verticalLayout_2 = QVBoxLayout(self.right_label_page)
        self.right_verticalLayout_2.setContentsMargins(3, 3, 3, 3)
        self.right_verticalLayout_2.setSpacing(2)

        # 上半部分
        self.right_function_region = QWidget()

        self.right_function_verticalLayout_21 = QVBoxLayout(self.right_function_region)
        self.right_function_verticalLayout_21.setContentsMargins(3, 3, 3, 3)
        self.right_function_verticalLayout_21.setSpacing(2)
        self.right_function_scrollArea_21 = QScrollArea(self.right_function_region)
        self.right_function_scrollArea_21.setFrameShape(QFrame.Panel)
        self.right_function_scrollArea_21.setFrameShadow(QFrame.Raised)
        self.right_function_scrollArea_21.setWidgetResizable(True)
        self.right_function_scrollArea_21.verticalScrollBar().setStyleSheet(scroll_str)
        self.right_function_scrollArea_21.horizontalScrollBar().setStyleSheet(scroll_str)
        self.right_function_scrollAreaWidgetContents_21 = QWidget()
        self.right_function_scrollAreaWidgetContents_21.setGeometry(QRect(0, 0, 348, 393))

        self.right_function_verticalLayout = QVBoxLayout(self.right_function_scrollAreaWidgetContents_21)

        self.right_tic = QWidget()
        self.right_tic_verticalLayout = QVBoxLayout(self.right_tic)
        self.right_tic_verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.right_tic_verticalLayout.setSpacing(10)
        self.right_config = QWidget()
        self.right_config_verticalLayout = QVBoxLayout(self.right_config)
        self.right_config_verticalLayout.setSpacing(2)

        self.ppm_region = QWidget()
        self.ppm_region_horizontalLayout = QHBoxLayout(self.ppm_region)
        self.ppm_region_horizontalLayout.setSpacing(4)
        self.ppm_region.setFixedWidth(200)
        self.ppm_info = QLabel()
        self.ppm_info.setText("ppm")
        self.ppm_info.setStyleSheet("QLabel{background:none}")
        self.ppm_info.setFixedWidth(36)
        self.edit_ppm = QLineEdit()
        self.edit_ppm.setStyleSheet("QLabel{background:none}")
        self.edit_ppm.setPlaceholderText(str(self.ppm))
        self.edit_ppm.textEdited.connect(self.ppm_text_changed)
        self.edit_ppm.setObjectName("font_gray")
        self.edit_ppm.setFixedWidth(50)
        self.right_anylyse = QPushButton("应用")
        self.right_anylyse.clicked.connect(self.apply_ppm)
        self.right_anylyse.setFixedWidth(58)
        self.right_anylyse.setStyleSheet(scan_btn_str)
        self.ppm_region_horizontalLayout.addWidget(self.ppm_info)
        self.ppm_region_horizontalLayout.addWidget(self.edit_ppm)
        self.ppm_region_horizontalLayout.addWidget(self.right_anylyse)

        self.right_config_verticalLayout.addWidget(self.right_tic)
        self.right_config_verticalLayout.addWidget(self.ppm_region)

        self.right_function_verticalLayout.addWidget(self.right_config)
        self.right_function_verticalLayout.addStretch()
        self.right_function_scrollArea_21.setWidget(self.right_function_scrollAreaWidgetContents_21)
        self.right_function_verticalLayout_21.addWidget(self.right_function_scrollArea_21)

        self.right_center = QWidget()
        self.right_center_horizontalLayout = QVBoxLayout(self.right_center)
        self.right_center_horizontalLayout.setContentsMargins(1, 1, 1, 1)
        self.right_center_horizontalLayout.setSpacing(2)
        self.right_center_horizontalLayout.setAlignment(Qt.AlignHCenter)
        # 下半部分
        self.right_struct_info = QWidget()

        self.right_verticalLayout_21 = QVBoxLayout(self.right_struct_info)
        self.right_verticalLayout_21.setContentsMargins(3, 3, 3, 3)
        self.right_verticalLayout_21.setSpacing(2)
        self.right_scrollArea_21 = QScrollArea(self.right_struct_info)
        self.right_scrollArea_21.setFrameShape(QFrame.Panel)
        self.right_scrollArea_21.setFrameShadow(QFrame.Raised)
        self.right_scrollArea_21.setWidgetResizable(True)
        self.right_scrollArea_21.horizontalScrollBar().setStyleSheet(scroll_str)
        self.right_scrollArea_21.verticalScrollBar().setStyleSheet(scroll_str)
        self.right_scrollAreaWidgetContents_21 = QWidget()
        self.right_scrollAreaWidgetContents_21.setGeometry(QRect(0, 0, 348, 393))

        self.right_verticalLayout_211 = QVBoxLayout(self.right_scrollAreaWidgetContents_21)
        self.right_verticalLayout_211.setSpacing(2)
        self.right_verticalLayout_211.setObjectName("right_verticalLayout_20")

        self.struct_title_region = QWidget()
        self.struct_title_region_horizontalLayout = QHBoxLayout(self.struct_title_region)

        self.struct_title_num = QLabel()
        self.struct_title_info = QLabel()
        self.struct_title_score = QLabel()
        self.struct_title_region_horizontalLayout.addWidget(self.struct_title_num)
        self.struct_title_region_horizontalLayout.addWidget(self.struct_title_info)
        self.struct_title_region_horizontalLayout.addWidget(self.struct_title_score)

        self.struct_title_info_2 = QLabel()
        self.struct_title_info_2.setObjectName("font8")

        self.right_comp_region = QWidget()
        self.right_label_horizontalLayout = QHBoxLayout(self.right_comp_region)
        self.right_label_horizontalLayout.setSpacing(2)

        self.right_struct_btn_region = QWidget()
        self.right_struct_btn_verticalLayout = QVBoxLayout(self.right_struct_btn_region)
        self.right_struct_region = QWidget()
        self.right_struct_region_verticalLayout = QVBoxLayout(self.right_struct_region)
        self.right_struct_score = QWidget()
        self.right_struct_score_verticalLayout = QVBoxLayout(self.right_struct_score)
        self.right_struct_btn_verticalLayout.setSpacing(2)
        self.right_struct_region_verticalLayout.setSpacing(2)
        self.right_struct_score_verticalLayout.setSpacing(2)
        self.right_label_horizontalLayout.addWidget(self.right_struct_btn_region)
        self.right_label_horizontalLayout.addWidget(self.right_struct_region)
        self.right_label_horizontalLayout.addWidget(self.right_struct_score)

        self.right_verticalLayout_211.addWidget(self.struct_title_region)
        self.right_verticalLayout_211.addWidget(self.struct_title_info_2)
        self.right_verticalLayout_211.addWidget(self.right_comp_region)

        self.right_verticalLayout_211.addStretch()

        self.right_scrollArea_21.setWidget(self.right_scrollAreaWidgetContents_21)
        self.right_verticalLayout_21.addWidget(self.right_scrollArea_21)

        sub_splitter = QSplitter(self)
        sub_splitter.setOrientation(Qt.Vertical)
        sub_splitter.addWidget(self.right_function_region)
        sub_splitter.addWidget(self.right_center)
        sub_splitter.addWidget(self.right_struct_info)  # 左侧FigureCanvas对象

        sub_splitter.setStretchFactor(0, 40)
        sub_splitter.setStretchFactor(1, 10)
        sub_splitter.setStretchFactor(2, 50)
        sub_splitter.setHandleWidth(1)
        self.right_verticalLayout_2.addWidget(sub_splitter)

        # 左侧面板
        self.right_label_page.setStyleSheet(qstr)
        self.right_label_page.setFixedWidth(240)

        # 中间图表面板
        self._fig = mpl.figure.Figure(figsize=(8, 5), dpi=72)  # 单位英寸
        self._fig.subplots_adjust(left=0.09, bottom=0.09, right=0.95, top=0.9, wspace=0.3, hspace=0.3)
        self.figCanvas = FigureCanvas(self._fig)  # 创建FigureCanvas对象，必须传递一个Figure对象
        self.figCanvas.setCursor(Qt.PointingHandCursor)
        self.figCanvas.setStyleSheet("QWidget{border:none}")

        self._drawInit()
        self._draw3D_init()

        # 设置图的工具栏
        self._initFigToolBar()

        # 设置左侧边栏和右侧图表面板
        self.fig_splitter = QSplitter(self)
        self.fig_splitter.setOrientation(Qt.Vertical)
        self.fig_splitter.addWidget(self.figCanvas)
        self.fig_splitter.addWidget(self.naviToolbar)
        self.fig_splitter.setStretchFactor(0, 92)
        self.fig_splitter.setStretchFactor(1, 8)

        self.splitter = QSplitter(self)
        self.splitter.setOrientation(Qt.Horizontal)
        self.splitter.addWidget(self.right_label_page)
        self.splitter.addWidget(self.fig_splitter)
        # if platform.system() != "Windows":
        #     leftWidth = int(660 * 100.0 / (QDesktopWidget().screenGeometry().width()))
        # else:
        #     leftWidth = int(840 * 100.0 / (QDesktopWidget().screenGeometry().width()))
        # self.splitter.setStretchFactor(0, leftWidth)
        # self.splitter.setStretchFactor(1, 100 - leftWidth)
        self.splitter.setHandleWidth(0)
        self.splitter.setStyleSheet("QSplitter::handle { background-color:white}")
        self.setCentralWidget(self.splitter)
        self.right_label_page.hide()
        self.naviToolbar.hide()

    def _initMenu(self):

        # 菜单栏的action
        openMzFileAct = QAction('Open .mzML', self)
        openMzFileAct.setStatusTip('Open .mzML')
        openMzFileAct.triggered.connect(self.openTICByMzML)

        openRawFileAct = QAction('Open .raw', self)
        openRawFileAct.setStatusTip('Open .raw')
        openRawFileAct.triggered.connect(self.openTICByRawFile)

        openDirAct = QAction('Open raw dir', self)
        openDirAct.setStatusTip('Open raw directory')
        openDirAct.triggered.connect(self.openTICByRawDir)

        showTICAct = QAction('show TIC', self)
        showTICAct.triggered.connect(self.showTIC)

        showTableAct = QAction('show table', self)
        showTableAct.triggered.connect(self._drawTable)

        downloadAct = QAction('download table', self)
        downloadAct.setStatusTip('download table')
        downloadAct.triggered.connect(self.downloadTabel)

        aboutAct = QAction('关于hepParser', self)
        aboutAct.setStatusTip('about hepParser')
        aboutAct.triggered.connect(self.aboutUs)

        exitMenuAct = QAction('Exit', self)
        exitMenuAct.setStatusTip('Exit application')
        exitMenuAct.triggered.connect(self.close)

        # 菜单栏
        menubar = self.menuBar()

        fileMenu = menubar.addMenu('File')  # File还是&File
        fileMenu.addAction(openMzFileAct)
        fileMenu.addAction(openRawFileAct)
        fileMenu.addAction(openDirAct)

        viewMenu = menubar.addMenu('View')
        viewMenu.addAction(showTICAct)
        viewMenu.addAction(showTableAct)

        toolMenu = menubar.addMenu('Tools')
        toolMenu.addAction(downloadAct)

        aboutMenu = menubar.addMenu('about')
        aboutMenu.addAction(aboutAct)

        exitMenu = menubar.addMenu('Exit')
        exitMenu.addAction(exitMenuAct)

    def _initToolBar(self):
        # 自定义的工具栏工具
        self.openDirAction = QAction(QIcon('icon/folder.svg'), "打开目录", self)
        self.openDirAction.setStatusTip('Open raw directory')
        self.openDirAction.triggered.connect(self.openTICByRawDir)
        self.openDirAction.setObjectName("blue")

        self.openFileAction = QAction(QIcon('icon/file.svg'), "打开文件", self)
        self.openFileAction.setStatusTip('Open .mzml')
        self.openFileAction.triggered.connect(self.openTICByMzML)

        self.showTICAction = QAction(QIcon('icon/T.svg'), 'TIC图', self)
        self.showTICAction.setStatusTip('TIC图')
        self.showTICAction.triggered.connect(self.showTIC)

        self.changeAction = QAction(QIcon('icon/report.svg'), '修改组成', self)
        self.changeAction.setStatusTip('修改组成')
        self.changeAction.triggered.connect(self.change_comp)

        self.tableAction = QAction(QIcon('icon/form.svg'), '显示表格', self)
        self.tableAction.setStatusTip('显示表格')
        self.tableAction.triggered.connect(self._drawTable)

        self.downloadAction = QAction(QIcon('icon/download.svg'), '下载组成', self)
        self.downloadAction.setStatusTip('下载组成')
        self.downloadAction.triggered.connect(self.downloadTabel)

        # 设置工具栏的语言
        self.tool_bar = QToolBar()
        self.tool_bar.setMovable(False)
        self.tool_bar.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)
        self.tool_bar.addAction(self.openDirAction)
        self.tool_bar.addAction(self.openFileAction)
        self.tool_bar.addAction(self.showTICAction)
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        spacer.setEnabled(False)
        self.tool_bar.addWidget(spacer)
        self.tool_bar.addAction(self.changeAction)
        self.tool_bar.addAction(self.tableAction)
        self.tool_bar.addAction(self.downloadAction)

        self.tool_bar.setStyleSheet(navi_str)
        self.tool_bar.setIconSize(QSize(24, 24))
        self.addToolBar(Qt.LeftToolBarArea, self.tool_bar)

    def _initFigToolBar(self):
        self.naviToolbar = NavigationToolbar(self.figCanvas, self)  # 创建工具栏
        #self.naviToolbar.setToolButtonStyle(Qt.ToolButtonTextUnderIcon)

        self.naviToolbar.setToolButtonStyle(Qt.ToolButtonIconOnly)
        actList = self.naviToolbar.actions()  # 关联的Action列表
        # actList[0].setText("复位")  # Home
        actList[0].setText("")  # Home
        actList[0].setToolTip("复位到原始视图")  # Reset original view
        actList[0].setIcon(QIcon('icon/home-run.svg'))

        actList[1].setText("")  # Back
        # actList[1].setText("回退")  # Back
        actList[1].setToolTip("回到前一视图")  # Back to previous view
        actList[1].setIcon(QIcon('icon/left-arrow.svg'))

        # actList[2].setText("前进")  # Forward
        actList[2].setText("")  # Forward
        actList[2].setToolTip("前进到下一视图")  # Forward to next view
        actList[2].setIcon(QIcon('icon/next.svg'))

        actList[4].setText("")  # Pan
        # actList[4].setText("拖拽")  # Pan
        actList[4].setToolTip("左键平移坐标轴，右键缩放坐标轴")  # Pan axes with left mouse, zoom with right
        actList[4].setIcon(QIcon('icon/move.svg'))

        # actList[5].setText("缩放")  # Zoom
        actList[5].setText("")  # Zoom
        actList[5].setToolTip("框选矩形框缩放")  # Zoom to rectangle
        actList[5].setIcon(QIcon('icon/zoom.svg'))

        actList[6].setText("")  # Subplots
        # actList[6].setText("子图")  # Subplots
        actList[6].setToolTip("设置子图")  # Configure subplots
        actList[6].setIcon(QIcon('icon/desktop.svg'))

        # actList[7].setText("定制")  # Customize
        actList[7].setText("")  # Customize
        actList[7].setToolTip("定制图表参数")  # Edit axis, curve and image parameters
        actList[7].setIcon(QIcon('icon/config.svg'))

        # actList[9].setText("保存图片")  # Save
        actList[9].setText("")  # Save
        actList[9].setToolTip("保存图表")  # Save the figure
        actList[9].setIcon(QIcon('icon/save.svg'))
        # 设置初始透明度:
        self.setOpacity()
        self.naviToolbar.insertAction(actList[0], actList[8])

        self.naviToolbar.removeAction(actList[3])
        self.naviToolbar.removeAction(actList[6])
        self.naviToolbar.removeAction(actList[7])
        self.naviToolbar.removeAction(actList[8])
        self.naviToolbar.removeAction(actList[10])
        self.naviToolbar.setIconSize(QSize(18, 18))

        self.naviToolbar.setStyleSheet(fig_navi_str)

    def _drawInit(self):
        img = Image.open('icon/background.png')
        ax1 = self._fig.add_subplot(1, 1, 1)
        ax1.imshow(img)
        ax1.axis("off")
        self._fig.canvas.draw_idle()

    def generate_left_side(self, flag='range'):
        self.right_label_page.show()

        tmp_region = QWidget()
        tmp_region.setFixedWidth(180)
        tmp_hor = QHBoxLayout(tmp_region)
        tmp_hor.setContentsMargins(1, 1, 1, 1)
        tmp_hor.setSpacing(8)

        if flag == 'range':
            tmp_btn_1 = QPushButton(str(self.start_scan) + "-" + str(self.end_scan))
        else:
            tmp_btn_1 = QPushButton(str(self.scan))

        tmp_btn_1.setStyleSheet(tmp_scan_label_str)
        tmp_btn_1.setFixedWidth(100)

        tmp_btn_2 = QPushButton(" ")
        tmp_btn_2.clicked.connect(self.hep_analyse)
        tmp_btn_2.setStyleSheet(scan_icon_str)
        tmp_btn_2.setFixedWidth(30)

        tmp_hor.addWidget(tmp_btn_1)
        tmp_hor.addWidget(tmp_btn_2)
        self.right_tic_verticalLayout.addWidget(tmp_region)

    def load_merged_peaks(self):
        self.peaks, self.maxIntensity, self.total_int, self.exp_isp, self.the_spectra, self.dict_list = get_filter_MZ(
            origin_MZ=self.spectrum,
            max_int=self.maxIntensity,
            fit_list=self.fit_list,
            the_HP=self.the_HP,
            ppm=self.ppm,
            bound_th=self.bound_th,
            bound_intensity=self.bound_intensity)
        for peak in self.peaks:
            self.peak_dict[peak[0]] = peak[1]
        self.origin_xs, self.origin_ys = format_peaks(self.peaks)

    def plot_origin_peaks(self, xs, ys):  # 获取原始峰
        self.naviToolbar.show()
        self.orgAx.plot(xs, ys, linewidth=1, c='gray', zorder=1)
        self.orgAx.set_ybound(lower=0, upper=1.1 * max(ys))
        self.orgAx.set_xlabel('m/z', fontsize=18)
        self.orgAx.set_ylabel('intensity', fontsize=18)
        self.orgAx.spines['top'].set_visible(False)
        self.orgAx.spines['right'].set_visible(False)
        self.orgAx.spines['left'].set_linewidth(1.5)
        self.orgAx.spines['bottom'].set_linewidth(1.5)
        self.orgAx.set_title("$MS^1$")  # 子图标题
        self._fig.canvas.draw_idle()
        self.opacity = [True, False, False, False]
        self.setOpacity()

    def _analyse_Composition(self):

        # 获取全部组成
        self.dataDlgProcess.close()
        self.all_right_struct, self.all_key_with_order = get_all_struct(self.label_info)
        self.chk_list = []
        self.change_comp()

    def _draw_composition(self):
        # 标注分析
        self.naviToolbar.show()
        self.figFlag[0] = 3
        self.load_label()
        self.mass_anotation = {}

        self._fig.clear()
        self.orgAx = self._fig.add_subplot(1, 1, 1)

        # 删除按钮
        for i in range(self.right_center_horizontalLayout.count()):
            self.right_center_horizontalLayout.itemAt(i).widget().deleteLater()
        for i in range(self.right_struct_btn_verticalLayout.count()):
            self.right_struct_btn_verticalLayout.itemAt(i).widget().deleteLater()
            self.right_struct_region_verticalLayout.itemAt(i).widget().deleteLater()
            self.right_struct_score_verticalLayout.itemAt(i).widget().deleteLater()

        # 重新添加按钮
        change_btn = QPushButton("启动单分子标注")
        change_info_region = QWidget()
        change_info_region_hori = QVBoxLayout(change_info_region)
        change_info_region_hori.setAlignment(Qt.AlignHCenter)
        change_info_region_hori.setContentsMargins(0, 0, 0, 0)
        change_info = QLabel()
        change_info.setText("点击后可选择单个分子对谱进行解释")
        change_info.setStyleSheet("QLabel{font-size:8px;}")
        change_info_region_hori.addWidget(change_info)
        change_btn.clicked.connect(self._label)
        change_btn.setFixedWidth(150)
        self.right_center_horizontalLayout.addWidget(change_btn)
        self.right_center_horizontalLayout.addWidget(change_info_region)
        self.struct_title_num.setText("&nbsp;序号<sup></sup>")
        self.struct_title_info.setText("&nbsp;&nbsp;&nbsp;分子组成<sup>a</sup>")
        self.struct_title_info_2.setText("&nbsp;&nbsp;&nbsp;&nbsp;<sup>a</sup>[HexA,GlcA,GlcN,Ac,SO3,Levoglucosan,Man]")
        self.struct_title_score.setText("<sup></sup>&nbsp;&nbsp;&nbsp;&nbsp;单分子解释度")
        for i in range(len(self.right_struct)):
            label_btn = QPushButton(str(i + 1))
            label_btn.setStyleSheet(
                label_btn_str[:-1] + "QPushButton{background-color:" + self.colors_str[self.right_struct[i][2]] + "}")

            label_btn.setFixedWidth(21)

            label_struct = QPushButton(self.right_struct[i][0])
            label_struct.setStyleSheet(label_struct_str)

            label_score = QPushButton(str(self.right_struct[i][1]))
            label_score.setStyleSheet(label_struct_str)

            label_btn.clicked.connect(lambda: self._labelFamilyPeak(self.sender().text()))
            label_struct.clicked.connect(lambda: self._labelFamilyPeak(self.struct_id[self.sender().text()]))

            self.right_struct_btn_verticalLayout.addWidget(label_btn)
            self.right_struct_region_verticalLayout.addWidget(label_struct)
            self.right_struct_score_verticalLayout.addWidget(label_score)

        for i in range(self.right_struct_btn_verticalLayout.count()):
            self.right_struct_btn_verticalLayout.itemAt(i).widget().setDisabled(True)
            self.right_struct_region_verticalLayout.itemAt(i).widget().setDisabled(True)
            self.right_struct_score_verticalLayout.itemAt(i).widget().setDisabled(True)
        self.massToStructs = {}
        self.massToStructsScore = {}
        self.massToStructsChks = {}
        for key in self.xy:
            label_list = self.mass_labels[key]
            a_list = []

            for j in range(len(label_list)):
                if label_list[j][4] != 1:
                    self.orgAx.annotate(s="+" + str(label_list[j][4] - 1), xy=key,
                                        xytext=(-self.x_space, +(10 + self.y_space * j + self.y_space * 0.5)),
                                        textcoords='offset pixels',
                                        color='black', ha='center', va='bottom',
                                        fontsize=8, fontweight='extra bold')

                a = self.orgAx.annotate(s=str(label_list[j][-2]), xy=key,
                                        xytext=(+0, +(10 + self.y_space * j)),
                                        textcoords='offset pixels',
                                        color='white', ha='center', va='bottom',
                                        fontsize=self.anno_text_size, fontweight='extra bold',
                                        bbox=dict(boxstyle='circle,pad=0.1',
                                                  fc=self.colors[label_list[j][-2]],
                                                  ec=self.colors[label_list[j][-2]],
                                                  lw=1))
                self.orgAx.annotate(s=str(label_list[j][2]) + "-", xy=key,
                                    xytext=(+self.x_space, +(10 + self.y_space * j + self.y_space * 0.5)),
                                    textcoords='offset pixels',
                                    color='black', ha='center', va='bottom',
                                    fontsize=8, fontweight='extra bold')
                a_list.append(a)
            self.mass_anotation[key] = a_list
        # 画原始结构
        artists = []
        for i in range(len(self.label_xs)):
            a = self.orgAx.plot(self.label_xs[i], self.label_ys[i],
                                label=self.mass_struct_tips[(self.label_xs[i][1], self.label_ys[i][1])],
                                linewidth=1, color='blue', picker='True', zorder=2)
            artists.append(a)
        mplcursors.cursor(pickables=self.orgAx, hover=True,
                          annotation_kwargs=dict(
                              bbox=dict(boxstyle="round,pad=.7",
                                        fc="#e3e3e3",
                                        alpha=.2,
                                        ec="#7b94b1"),
                              arrowprops=dict(
                                  arrowstyle="->",
                                  connectionstyle="arc3",
                                  shrinkB=3,
                                  ec="#7b94b1",
                              ))).connect("add", lambda sel: sel.annotation.set_text(sel.artist.get_label()))
        self.orgAx.plot(self.origin_xs, self.origin_ys, linewidth=1, c='gray', zorder=1)
        self.figFlag[0] = 3

        self.orgAx.set_ybound(lower=0)
        self.orgAx.set_xlabel('m/z', fontsize=18)
        self.orgAx.set_ylabel('intensity', fontsize=18)
        self.orgAx.spines['top'].set_visible(False)
        self.orgAx.spines['right'].set_visible(False)
        self.orgAx.spines['left'].set_linewidth(1.5)
        self.orgAx.spines['bottom'].set_linewidth(1.5)
        self.orgAx.set_title("$MS^1$")  # 子图标题
        self._fig.canvas.draw_idle()

        self.opacity = [True, True, True, True]
        self.setOpacity()
        self.tableAction.setText("显示表格")
        self.tableAction.setIcon(QIcon('icon/form.svg'))
        self.tableAction.triggered.disconnect()
        self.tableAction.triggered.connect(self._drawTable)
        self.right_anylyse.setStyleSheet(scan_btn_str)

    def _drawTable(self):
        self.naviToolbar.hide()
        if self.figFlag[0] < 3:
            QMessageBox.information(self, "Message", "请先加载数据，然后选谱分析")
        else:
            self._fig.clear()
            ax2 = self._fig.add_subplot(1, 1, 1)
            rowLabel = ['质荷比', '电荷', '第几个同位素峰', '分子组成', '脱落基团']

            table = ax2.table(cellText=self.labels[0:20], colLabels=rowLabel, rowLoc='center',
                              loc='center', cellLoc='center', fontsize=60, edges='open')
            for (row, col), cell in table.get_celld().items():
                if row == 0:
                    cell.visible_edges = "BT"
                if row == 20:
                    cell.visible_edges = 'B'
            # ax2.annotate('目前只显示前18行数据，更多数据请点击"下载组成"下载后查看', xy=(0.6, 0), color='black', va='bottom', fontsize=10)
            ax2.set_title(
                '目前只显示前20行数据，更多数据请点击"下载组成"下载后查看\n分子组成:[$HexA,GlcA,GlcN,Ac,SO3,Levoglucosan,Man$],基团脱落:[$HSO_3$, $NH_2$, $NH_{2}SO_3$, $COOH$]')
            ax2.axis("off")
            table.scale(1, 2.3)

        self._fig.canvas.draw_idle()
        self.tableAction.setText("返回谱图")
        self.tableAction.setIcon(QIcon('icon/back.svg'))
        self.tableAction.triggered.disconnect()
        self.tableAction.triggered.connect(self._draw_composition)

    def _label(self):
        self.naviToolbar.show()
        if self.figFlag[0] < 3:
            QMessageBox.information(self, "Message", "请先加载数据，然后选谱分析")
        else:
            self.figFlag[0] = 4
            self._labelFamilyPeak(1)
            # 重新添加按钮
            for i in range(self.right_center_horizontalLayout.count()):
                self.right_center_horizontalLayout.itemAt(i).widget().deleteLater()
            self.struct_title_num.setText("&nbsp;序号<sup></sup>")
            self.struct_title_info.setText("&nbsp;&nbsp;&nbsp;分子组成<sup>a</sup>")
            self.struct_title_info_2.setText(
                "&nbsp;&nbsp;&nbsp;&nbsp;<sup>a</sup>[HexA,GlcA,GlcN,Ac,SO3,Levoglucosan,Man]")
            self.struct_title_score.setText("<sup></sup>&nbsp;&nbsp;&nbsp;&nbsp;单分子解释度")

            change_btn = QPushButton("关闭单分子标注")  # 启动单分子标注
            change_info_region = QWidget()
            change_info_region_hori = QVBoxLayout(change_info_region)
            change_info_region_hori.setAlignment(Qt.AlignHCenter)
            change_info_region_hori.setContentsMargins(0, 0, 0, 0)
            change_info = QLabel()
            change_info.setText("返回谱图")
            change_info.setStyleSheet("QLabel{font-size:8px;}")
            change_info_region_hori.addWidget(change_info)

            change_btn.setFixedWidth(150)
            change_btn.clicked.connect(self._draw_composition)
            self.right_center_horizontalLayout.addWidget(change_btn)
            self.right_center_horizontalLayout.addWidget(change_info_region)
            for i in range(self.right_struct_btn_verticalLayout.count()):
                self.right_struct_btn_verticalLayout.itemAt(i).widget().setDisabled(False)
                self.right_struct_region_verticalLayout.itemAt(i).widget().setDisabled(False)
                self.right_struct_score_verticalLayout.itemAt(i).widget().setDisabled(False)

    def _labelFamilyPeak(self, index):
        for i in range(self.right_struct_btn_verticalLayout.count()):
            self.right_struct_btn_verticalLayout.itemAt(i).widget().setStyleSheet(
                label_btn_str[:-1] + "QPushButton{background-color:" + self.colors_str[self.right_struct[i][2]] + "}")
            self.right_struct_btn_verticalLayout.itemAt(i).widget().setFixedWidth(21)
            self.right_struct_region_verticalLayout.itemAt(i).widget().setStyleSheet(label_struct_str)
            self.right_struct_score_verticalLayout.itemAt(i).widget().setStyleSheet(label_struct_str)

        index = int(index) - 1
        self.right_struct_btn_verticalLayout.itemAt(index).widget().setStyleSheet(push_label_btn_str)
        self.right_struct_region_verticalLayout.itemAt(index).widget().setStyleSheet(push_label_struct_str)
        self.right_struct_score_verticalLayout.itemAt(index).widget().setStyleSheet(push_label_struct_str)

        self._fig.clear()
        self.ax1 = self._fig.add_subplot(1, 1, 1)

        self.ax1.plot(self.origin_xs, self.origin_ys, linewidth=1, c='gray')  # 用折线图画柱体

        # 标注出衍生峰和基团损失

        mz_list, inten_list, comp_list, lose_list, z_list, score_list = self.mass_family[self.key_with_order[index]]
        xs, ys = format_peaks_2_dim(mzs=mz_list, intens=inten_list)
        self.ax1.plot(xs, ys, linewidth=1, c='blue', label=self.right_struct[index][0])

        for i in range(len(mz_list)):
            lose_info = format_loss(lose_list[i], self.lose_ion)
            if lose_info == "$$":
                lose_info = format_comp(comp_list[i], z_list[i])
            self.ax1.annotate(s=lose_info, xy=(mz_list[i], inten_list[i]), xytext=(+0, +20),
                              color='blue', textcoords='offset pixels', rotation=90, ha='center', va='bottom',
                              fontsize=11)

        self.ax1.set_ybound(lower=0)
        self.ax1.set_xlabel('m/z', fontsize=18)
        self.ax1.set_ylabel('intensity', fontsize=18)
        self.ax1.spines['top'].set_visible(False)
        self.ax1.spines['right'].set_visible(False)
        self.ax1.spines['left'].set_linewidth(1.5)
        self.ax1.spines['bottom'].set_linewidth(1.5)
        self.ax1.set_title(("$MS^1$"))
        self.ax1.legend()
        self._fig.canvas.draw_idle()

        self.opacity = [True, False, False, True]
        self.setOpacity()

    def load_label(self):  # 转换标注数据

        # 3.标注匹配上结构，红色
        print("family:", len(self.label_info[3]))
        #  self.labels = format_labels = self.label_info[0]

        self.xy, self.mass_labels, self.mass_struct_tips, self.right_struct, self.struct_id, self.label_message, self.labels \
            = get_labels(self.label_info, self.peak_dict, self.lose_ion, self.key_with_order)

        self.label_xs, self.label_ys = format_peaks_alone(self.xy)

        # 4.衍生峰里的
        self.mass_family = get_family(self.label_info, self.peak_dict)

        # 颜色生成器
        colors = mpl.cm.get_cmap("YlOrRd")  # 'viridis', 'inferno', 'plasma', 'magma'
        digit = list(map(str, range(10))) + list("ABCDEF")
        self.colors = colors(np.linspace(0, 1, len(self.right_struct) + 2))
        self.colors_str = []
        for color in self.colors:
            string = '#'
            for i in range(0, 3):
                a1 = int(color[i] * 255) // 16
                a2 = int(color[i] * 255) % 16
                string += digit[a1] + digit[a2]
            self.colors_str.append(string)
        self.colors = self.colors.tolist()
        self.colors.reverse()
        self.colors_str.reverse()

    def set_xmlInfo(self, xmlInfo):
        self.xmlInfo = xmlInfo

    def updateConvertProcessBar(self, preMess):
        curr, total = preMess.split('/')
        QApplication.processEvents()  # 实时刷新界面
        self.dlgProgress.setValue(int(float(curr) * 100.0 / float(total)))
        QApplication.processEvents()  # 实时刷新界面
        self.dlgProgress.setLabelText("正在转换质谱 " + preMess)
        QApplication.processEvents()  # 实时刷新界面

    def updateDataProcessBar(self, ID):
        QApplication.processEvents()  # 实时刷新界面
        self.dataDlgProcess.setValue(int(ID * 100.0 / self.total_comp))
        QApplication.processEvents()  # 实时刷新界面
        self.dataDlgProcess.setLabelText("已完成 " + str(ID) + "/" + str(self.total_comp) + "       ")
        QApplication.processEvents()  # 实时刷新界面

    def update_prpcess_step_bar(self, step):
        QApplication.processEvents()  # 实时刷新界面
        if step == 0:
            self.data_process_bar.setLabelText("正 在 转 换 质 谱...")
        elif step == 1:
            self.data_process_bar.setLabelText("正 在 整 合 数 据...")
        elif step == 2:
            self.data_process_bar.setLabelText("正 在 平 滑 数 据...")
        QApplication.processEvents()  # 实时刷新界面
        if step == 3:
            self.data_process_bar.close()

    def infoParseProcess(self, data_info):
        self.match_result, self.label_info, self.candidate_max_num, self.candi_score = data_info

    def mxmlParseProcess(self):
        self.mzThread = mzMLWorker(xmlFileName=self.mzmlFileName)
        self.mzThread.sinXml.connect(self.set_xmlInfo)
        self.data_process_bar = DPstepBar()
        self.data_process_bar.show()
        self.mzThread.step.connect(self.update_prpcess_step_bar)
        self.mzThread.finished.connect(self._draw3D)
        self.mzThread.start()

    def getConvertProcess(self):  # msconvert 转换文件的进度条
        self.dlgProgress = QProgressDialog("正在转换质谱格式...", "", 1, 100, self)
        self.dlgProgress_hori = QVBoxLayout(self.dlgProgress)

        self.dlgProgress.setWindowTitle("转换质谱格式")
        self.dlgProgress.setWindowModality(Qt.WindowModal)  # 模态对话框
        self.dlgProgress.setCancelButton(None)  # 隐藏取消按钮
        self.dlgProgress.setAutoReset(True)  # calls reset() as soon as value() equals maximum()
        self.dlgProgress.setAutoClose(True)  # whether the dialog gets hidden by reset()
        self.dlgProgress.setMinimumDuration(0)
        self.dlgProgress.setFixedSize(280, 110)
        self.dlgProgress.setStyleSheet(pbar_str)
        self.dlgProgress.open()

        self.thread = ConvertWorker(curPath=self.curPath, rawFileName=self.rawFileName)
        self.thread.sinOut.connect(self.updateConvertProcessBar)
        self.thread.finished.connect(self.mxmlParseProcess)
        self.thread.start()

    def hep_analyse(self):
        print("debug:search")
        self.figFlag[0] = 3
        self.dataDlgProcess = QAnalyseBar()
        self.dataDlgProcess.show()
        self.dataThread = DataWorker(exp_isp=self.exp_isp, max_int=self.maxIntensity,
                                     total_int=self.total_int, the_spectra=self.the_spectra,
                                     dict_list=self.dict_list, the_HP=self.the_HP,
                                     ppm=self.ppm, bound_th=self.bound_th,
                                     bound_intensity=self.bound_intensity)
        self.dataThread.sinID.connect(self.updateDataProcessBar)
        self.dataThread.sinDataInfo.connect(self.infoParseProcess)
        self.dataThread.finished.connect(self._analyse_Composition)
        self.dataThread.start()

    def change_comp(self):
        if self.figFlag[0] < 2:
            QMessageBox.information(self, "Message", "请先加载数据，然后选谱分析")
        elif self.figFlag[0] == 2:
            QMessageBox.information(self, "Message", "请先解析谱，然后再进行混合物分析")
        else:
            if len(self.chk_list) == 0:
                better_index, better_score = get_first_max_num(self.candi_score)
                for i in range(0, len(self.candi_score)):
                    if i <= better_index:
                        self.chk_list.append(True)
                    else:
                        self.chk_list.append(False)
            dlgTableSize = QmyDialogSize(comps=self.all_right_struct, chk_list=self.chk_list,
                                         candi_score=self.candi_score)
            ret = dlgTableSize.exec()  # 模态方式运行对话框
            if (ret == QDialog.Accepted):
                self.chk_list = dlgTableSize.getCheckList()

            new_key_with_order = []
            for i in range(len(self.chk_list)):
                if self.chk_list[i]:
                    new_key_with_order.append(self.all_key_with_order[i])
            self.key_with_order = new_key_with_order
            self._draw_composition()

    def ppm_text_changed(self):
        self.right_anylyse.setStyleSheet(qstr)

    def apply_ppm(self):
        self.figFlag[0] = 3

        self.ppm = int(self.edit_ppm.text())
        # self.load_merged_peaks()
        self.dataDlgProcess = QAnalyseBar()
        self.dataDlgProcess.show()

        self.dataThread = DataWorker(exp_isp=self.exp_isp, max_int=self.maxIntensity, total_int=self.total_int,
                                     the_spectra=self.the_spectra, dict_list=self.dict_list,
                                     the_HP=self.the_HP, ppm=self.ppm, bound_th=self.bound_th,
                                     bound_intensity=self.bound_intensity)
        self.dataThread.sinID.connect(self.updateDataProcessBar)
        self.dataThread.sinDataInfo.connect(self.infoParseProcess)
        self.dataThread.finished.connect(self._analyse_Composition)
        self.dataThread.start()

    def openTICByMzML(self):
        mzmlFileName, fileType = QFileDialog.getOpenFileName(self, "打开文件", "", "*.mzML;;*.mzXML;;All Files(*)")
        if mzmlFileName != '':
            self.mzmlFileName = mzmlFileName
            self.mxmlParseProcess()

    def openTICByRawFile(self):
        rawFileName, fileType = QFileDialog.getOpenFileName(self, "打开文件", "", "*.raw;;*.RAW;;All Files(*)")

        if rawFileName != '':
            self.rawFileName = rawFileName
            self.curPath = os.path.abspath('.').replace("\\", "/") + '/mzML/'
            self.mzmlFileName = 'mzML/' + rawFileName.strip(".raw").split('/')[-1] + ".mzML"
            self.getConvertProcess()

    def openTICByRawDir(self):
        selectDir = QFileDialog.getExistingDirectory(self, "选择raw文件所在的目录", "", QFileDialog.ShowDirsOnly)

        if selectDir != '':
            self.rawFileName = selectDir
            self.curPath = os.path.abspath('.').replace("\\", "/") + '/mzML/'
            self.mzmlFileName = 'mzML/' + selectDir.split('/')[-1] + ".mzML"
            self.getConvertProcess()

    def showTIC(self):
        if self.figFlag[0] == 0:
            QMessageBox.information(self, "Message", "数据未加载,请先打开一个raw文件或目录")
        else:
            self._draw3D()

    def downloadTabel(self):
        if self.figFlag[0] < 3:
            QMessageBox.information(self, "Message", "请先加载数据，然后选谱分析")
        else:
            selectedDir, filtUsed = QFileDialog.getSaveFileName(self, "下载组成", 'data/formula.xls',
                                                                "*.xls;;All Files(*)")
            if selectedDir != '':
                self.save_xls(selectedDir)
                # self.df.to_csv(selectedDir)

    def save_xls(self, path):
        self.ws.write(0, 0, '质荷比')
        self.ws.write(0, 1, '电荷')
        self.ws.write(0, 2, '第几个同位素峰')
        self.ws.write(0, 3, '分子组成')
        self.ws.write(0, 4, '脱落基团')

        i = 1
        for key in self.xy:
            label_list = self.mass_labels[key]
            for j in range(len(label_list)):
                self.ws.write(i, 0, str(key[0]))
                self.ws.write(i, 1, str(label_list[j][2]))
                self.ws.write(i, 2, str(label_list[j][4]))
                self.ws.write(i, 3, str(label_list[j][0]))
                self.ws.write(i, 4, str(label_list[j][1]))
                i += 1

        # 保存excel文件
        self.wb.save(path)

    def setOpacity(self):  # 这是设置菜单栏4个按钮的背景色
        self.showTICAction.setEnabled(self.opacity[0])
        self.changeAction.setEnabled(self.opacity[1])
        self.tableAction.setEnabled(self.opacity[2])
        self.downloadAction.setEnabled(self.opacity[3])
        self.changeAction.setVisible(self.opacity[1])
        self.tableAction.setVisible(self.opacity[2])
        self.downloadAction.setVisible(self.opacity[3])

    def aboutUs(self):
        if platform.system() != "Windows":
            QMessageBox.about(self, "关于hepParser",
                              "     版本1.0.0(Beta版)      \n      Copyright © ICT       \n烟台中科生信智能科技中心\n中国科学院计算技术研究所\n                东诚药业")
        else:
            QMessageBox.about(self, "关于hepParser",
                              "     版本1.0.0(Beta版)      \n     Copyright © ICT       \n烟台中科生信智能科技中心\n中国科学院计算技术研究所\n            东诚药业")

    def closeEvent(self, event):
        reply = QMessageBox.question(self, 'Message',
                                     "Are you sure to quit?", QMessageBox.No | QMessageBox.Yes
                                     , QMessageBox.No)

        if reply == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

    # 画三维图
    def _draw3D_init(self):
        self.cursor = EventFactory(self._fig, self.figFlag)
        self._fig.canvas.mpl_connect('motion_notify_event', self.cursor.mouse_move)
        self._fig.canvas.mpl_connect('button_press_event', self.cursor.mouse_press)
        self._fig.canvas.mpl_connect('button_release_event', self.cursor.mouse_release)
        self._fig.canvas.mpl_connect('scroll_event', self.cursor.scroll)

    def _draw3D(self):
        self.right_label_page.hide()
        self.naviToolbar.hide()
        self.figFlag[0] = 1
        self.ppm = 20
        self._fig.clear()
        ax = self._fig.add_subplot(111, projection='3d')
        ax.disable_mouse_rotation()

        alldata = self.xmlInfo[-1]
        self.ori_mass_range = [alldata[0][0], alldata[0][-1]]
        self.ori_scan_range = [alldata[1][0], alldata[1][-1]]
        self.cursor.init(ax, alldata, self.ori_mass_range, self.ori_scan_range)
        for i in range(len(self.opacity)):
            self.opacity[i] = False
        self.setOpacity()
        for i in range(0, self.right_tic_verticalLayout.count()):
            self.right_tic_verticalLayout.itemAt(i).widget().deleteLater()

    # 创建右键菜单
    def _init_right_pick_menu(self):
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.show_context_menu)

        self.contextMenu = QMenu(self)
        self.merge_action = self.contextMenu.addAction(u'合并指定范围的质谱')
        self.singe_scan_action = self.contextMenu.addAction(u'选择单个时间点的质谱')

        self.cancel_action = self.contextMenu.addAction(u'取消')

        self.merge_action.triggered.connect(self.rightPick)
        self.singe_scan_action.triggered.connect(self.rightPick)
        self.cancel_action.triggered.connect(self.rightPickCancel)

    def show_context_menu(self, pos):
        '''''
        右键点击时调用的函数，菜单显示前，将它移动到鼠标点击的位置
        '''
        if self.figFlag[0]!=1:
            return
        self.contextMenu.move(self.pos() + pos)
        self.contextMenu.show()

    def rightPick(self):
        if self.figFlag[0] == 1:
            self.scan_range = self.cursor.get_current_scan_range()
            self.scan = self.cursor.get_current_scan()
            if self.scan_range is not None:
                self.start_scan, self.end_scan = self.scan_range

                # 生成左侧边栏
                self.generate_left_side(flag='range')

                # 加载未合并的峰
                self.spectrum, self.maxIntensity = get_merged_peaks(self.mzmlFileName, self.start_scan, self.end_scan)
                save_file(self.spectrum, "origin.mgf")

                # 画出原始图
                self.figFlag[0] = 2
                self._fig.clear()
                self.orgAx = self._fig.add_subplot(1, 1, 1)
                xs, ys = format_peaks(self.spectrum)
                save_file(self.peaks, "test.mgf")
                self.plot_origin_peaks(xs, ys)
                self.load_merged_peaks()
            elif self.scan is not None:
                self.spectrum, self.maxIntensity = get_peaks_by_id(self.mzmlFileName, self.scan)
                save_file(self.spectrum, "938.mgf")
                # 生成左侧边栏
                self.generate_left_side(flag='single')

                save_file(self.spectrum, "938.mgf")
                # 画出原始图
                self.figFlag[0] = 2
                self._fig.clear()
                self.orgAx = self._fig.add_subplot(1, 1, 1)
                xs, ys = format_peaks(self.spectrum)
                self.plot_origin_peaks(xs, ys)
                self.load_merged_peaks()
            else:
                QMessageBox.information(self, "Message", "未选择数据,请先选择一个scan或一段scan对应的信息")

    def rightPickCancel(self):
        pass


if __name__ == "__main__":
    # 每个PyQt5应用都必须创建一个应用对象。sys.argv是一组命令行参数的列表
    app = QApplication(sys.argv)
    splash = QSplashScreen(QPixmap("icon/loading.png"))
    splash.showMessage("", Qt.AlignHCenter | Qt.AlignBottom, Qt.black)
    splash.show()  # 显示启动界面
    qApp.processEvents()

    win = MainWindow()
    time.sleep(2)
    win.show()
    splash.finish(win)
    app.exit(app.exec_())
    dataClosed()