from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from spectra_tools import get_first_max_num

qstr = """
        QPushButton {
            border-radius: 4px;
            padding: 5px 12px;
            background-color: #3265f6;
            color: white;
        }
        QPushButton:hover {
            background-color: #3c7bf7;
        }
        QPushButton:pressed {
            background-color: #3c7bf7;
        }
    """

cancelStr = """
        QPushButton {
            border-radius: 4px;
            padding: 5px 12px;
            background-color: #dddddd;
            color: black;
        }
        QPushButton:hover {
            background-color: #cccccc;
        }
        QPushButton:pressed {
            background-color: #cccccc;
        }
    """

pbar_str = """
            QProgressBar {
                border:none;
                color:white;
                text-align:center;
                background-color:#e2e2e2;
                border-radius: 5px;
                height:20px;
                width:230px;
            }

            QProgressBar::chunk{
                border-radius:5px;
                background-color:#3265f6;
            }
            
            """


class Ui_QWDialogSize(object):

    def __init__(self, comps, chk_list, candi_score):
        self.comps = comps
        self.chk_list = chk_list
        self.candi_score = candi_score

    def setupUi(self, QWDialogSize):
        _translate = QtCore.QCoreApplication.translate

        QWDialogSize.setObjectName("QWDialogSize")
        QWDialogSize.setWindowTitle(_translate("QWDialogSize", "组成列表"))
        QWDialogSize.setWindowModality(QtCore.Qt.NonModal)
        QWDialogSize.resize(600, 400)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(QWDialogSize.sizePolicy().hasHeightForWidth())
        QWDialogSize.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setPointSize(10)
        QWDialogSize.setFont(font)
        QWDialogSize.setSizeGripEnabled(False)
        QWDialogSize.setModal(False)

        self.horizontalLayout = QtWidgets.QHBoxLayout(QWDialogSize)

        self.right_struct_info = QGroupBox(QWDialogSize)
        self.right_verticalLayout_21 = QVBoxLayout(self.right_struct_info)
        self.right_verticalLayout_21.setContentsMargins(3, 3, 3, 3)
        self.right_verticalLayout_21.setSpacing(10)
        self.right_scrollArea_21 = QScrollArea(self.right_struct_info)
        self.right_scrollArea_21.setFrameShape(QFrame.Panel)
        self.right_scrollArea_21.setFrameShadow(QFrame.Raised)
        self.right_scrollArea_21.setWidgetResizable(True)
        self.right_scrollAreaWidgetContents_21 = QWidget()
        self.right_scrollAreaWidgetContents_21.setGeometry(QRect(0, 0, 348, 393))

        self.right_verticalLayout_211 = QVBoxLayout(self.right_scrollAreaWidgetContents_21)
        self.right_verticalLayout_211.setContentsMargins(11, 11, 11, 11)
        self.right_verticalLayout_211.setSpacing(1)
        self.right_verticalLayout_211.setObjectName("right_verticalLayout_20")
        self.right_struct_title = QWidget()
        self.right_title_orizontalLayout = QHBoxLayout(self.right_struct_title)

        self.comp_label = QLabel(self.right_struct_info)
        self.comp_label.setText(_translate("MainWindow", '分子组成'))
        self.score_label = QLabel(self.right_struct_info)
        self.score_label.setText(_translate("MainWindow", '单分子解释度'))
        self.right_title_orizontalLayout.addWidget(self.comp_label)
        self.right_title_orizontalLayout.addWidget(self.score_label)

        self.right_comp_region = QWidget()
        self.right_label_horizontalLayout = QHBoxLayout(self.right_comp_region)
        self.right_label_horizontalLayout.setSpacing(0)

        self.right_struct_btn_region = QWidget()
        self.right_struct_btn_verticalLayout = QVBoxLayout(self.right_struct_btn_region)
        self.right_struct_score = QWidget()
        self.right_struct_score_verticalLayout = QVBoxLayout(self.right_struct_score)
        self.right_struct_btn_verticalLayout.setSpacing(1)
        self.right_struct_score_verticalLayout.setSpacing(1)
        self.right_label_horizontalLayout.addWidget(self.right_struct_btn_region)
        self.right_label_horizontalLayout.addWidget(self.right_struct_score)

        self.right_verticalLayout_211.addWidget(self.right_struct_title)
        self.right_verticalLayout_211.addWidget(self.right_comp_region)

        self.right_scrollArea_21.setWidget(self.right_scrollAreaWidgetContents_21)
        self.right_verticalLayout_21.addWidget(self.right_scrollArea_21)

        self.chkBox_list = []
        self.label_score_list = []
        for i in range(0, len(self.comps)):
            self.chkBox_list.append(QCheckBox(self.comps[i][0]))
            self.chkBox_list[i].setChecked(self.chk_list[i])

            self.label_score_list.append(QLabel())
            self.label_score_list[i].setText(str(self.comps[i][1]))

            self.right_struct_btn_verticalLayout.addWidget(self.chkBox_list[i])
            self.right_struct_score_verticalLayout.addWidget(self.label_score_list[i])

        self.horizontalLayout.addWidget(self.right_struct_info)

        self.frame = QFrame(QWDialogSize)
        self.frame.setMaximumSize(QtCore.QSize(90, 16777215))
        self.frame.setFrameShape(QFrame.StyledPanel)
        self.frame.setFrameShadow(QFrame.Raised)
        self.frame.setObjectName("frame")

        self.verticalLayout = QVBoxLayout(self.frame)
        self.verticalLayout.setObjectName("verticalLayout")

        self.selectAll_chk = QCheckBox("全选")
        self.selectAll_chk.setChecked(False)
        self.selectAll_chk.stateChanged.connect(self.select_all_compostion)
        self.verticalLayout.addWidget(self.selectAll_chk)

        self.select_addvise_chk = QCheckBox("推荐")
        self.select_addvise_chk.setChecked(True)
        self.select_addvise_chk.stateChanged.connect(self.select_advise)
        self.verticalLayout.addWidget(self.select_addvise_chk)

        self.btnOK = QPushButton(self.frame)
        self.btnOK.setObjectName("btnOK")
        self.btnOK.setText(_translate("QWDialogSize", "标注"))
        self.btnOK.setStyleSheet(qstr)
        self.verticalLayout.addWidget(self.btnOK)

        self.btnCancel = QPushButton(self.frame)
        self.btnCancel.setObjectName("btnCancel")
        self.btnCancel.setText(_translate("QWDialogSize", "返回"))
        self.btnCancel.setStyleSheet(cancelStr)
        self.verticalLayout.addWidget(self.btnCancel)

        self.horizontalLayout.addWidget(self.frame)

        self.btnOK.clicked.connect(QWDialogSize.accept)
        self.btnCancel.clicked.connect(QWDialogSize.reject)
        QtCore.QMetaObject.connectSlotsByName(QWDialogSize)

    def select_all_compostion(self):
        if self.selectAll_chk.isChecked():
            for i in range(0, len(self.candi_score)):
                self.chk_list[i] = True
                self.chkBox_list[i].setChecked(True)
        else:
            for i in range(0, len(self.candi_score)):
                self.chk_list[i] = False
                self.chkBox_list[i].setChecked(False)

    def select_advise(self):
        better_index, better_score = get_first_max_num(self.candi_score)
        if self.select_addvise_chk.isChecked():
            for i in range(0, len(self.candi_score)):
                if i <= better_index:
                    self.chk_list[i] = True
                    self.chkBox_list[i].setChecked(True)
                else:
                    self.chk_list[i] = False
                    self.chkBox_list[i].setChecked(False)
        else:
            for i in range(0, len(self.candi_score)):
                self.chk_list[i] = False
                self.chkBox_list[i].setChecked(False)


class QmyDialogSize(QDialog):
    def __init__(self, comps=None, chk_list=None, candi_score=None, parent=None):
        super().__init__(parent)  # 调用父类构造函数，创建窗体
        self.comps = comps
        self.chk_list = chk_list
        self.candi_score = candi_score
        self.ui = Ui_QWDialogSize(comps=self.comps, chk_list=self.chk_list, candi_score=self.candi_score)  # 创建UI对象
        self.ui.setupUi(self)  # 构造UI界面

        self.setWindowFlags(QtCore.Qt.MSWindowsFixedSizeDialogHint)

    def __del__(self):  # 析构函数
        print("QmyDialogSize 对象被删除了")

    def getCheckList(self):
        chk_list = []
        for i in range(0, len(self.comps)):
            chk_list.append(self.ui.chkBox_list[i].isChecked())
        if True not in chk_list:
            better_index, better_score = get_first_max_num(self.candi_score)
            if self.select_addvise_chk.isChecked():
                for i in range(0, len(self.candi_score)):
                    if i <= better_index:
                        self.chk_list[i] = True
        return chk_list


class QAnalyseBar(QDialog):

    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        self.setWindowTitle("肝素组成分析中")
        self.setWindowModality(Qt.NonModal)
        screen = QDesktopWidget().screenGeometry()
        windowSize = self.geometry()
        self.resize(260, 140)
        self.setStyleSheet("QWidget{background-color:white}")
        self.center = QWidget(self)
        self.center.setStyleSheet("QWidget{background-color:white}")
        self.verticalLayout = QVBoxLayout(self.center)
        self.horizontalLayout = QHBoxLayout()
        self.horizontalLayout.setContentsMargins(3, 3, 3, 10)
        self.horizontalLayout.setSpacing(0)

        self.gif_label = QLabel()
        self.precess_info = QLabel("即将开始组成分析")
        self.gif = QMovie('./icon/hep.gif')
        self.gif.setScaledSize(QSize(72, 60))

        self.gif_label.setMovie(self.gif)

        self.gif.start()

        self.pbar = QProgressBar(self)
        self.pbar.setStyleSheet(pbar_str)
        self.horizontalLayout.addWidget(self.gif_label)
        self.horizontalLayout.addWidget(self.precess_info)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.verticalLayout.addWidget(self.pbar)
        self.setModal(True)

    def setValue(self, value):
        self.pbar.setValue(value)

    def setLabelText(self, text):
        self.precess_info.setText(text)
