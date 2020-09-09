qstr = """
        QWidget {
            background:#fefefe;
            font-size:16px;
            color:#3265f6}
        QPushButton {
            border-radius: 4px;
            padding: 5px 12px;
            background-color: #3265f6;
            color: white;
        }
        QPushButton:hover {
            background-color: #538bf7;
        }
        QPushButton:pressed {
            background-color: #538bf7;
        }
        QWidget#left26{
            margin-left:20px;
        }
        QLineEdit{
            background-color: white;
            border:1px solid #dddddd;
            border-radius: 2px;
            padding: 6px 6px;
        }
        QMessageBox {
            font-size: 14px;
            padding: 10px;
        }
        QLabel#font8{
            font-size: 12px;
        }
        QLineEdit#font_gray{
            color:gray;
        }
        """
# 滚动条样式
scroll_str = """
            QScrollBar:vertical{
                width:6px;
                background:#dddddd;
                margin:0px,0px,0px,0px;
            }
            QScrollBar:horizontal{
                height:6px;
                background:#dddddd;
                margin:0px,0px,0px,0px;
            }
            QScrollBar::handle:vertical
            {
                width:6px;
                background:#bbbbbb; 
                border-radius:3px;   
                min-height:20;
            }
            QScrollBar::handle:horizontal
            {
                height:6px;
                background:#bbbbbb; 
                border-radius:3px;   
                min-width:20;
            }
            QScrollBar::handle:vertical:hover
            {
                width:6px;
                background:#aaaaaa;   
                border-radius:3px;
                min-height:20;
            }
            QScrollBar::handle:horizontal:hover
            {
                height:6px;
                background:#aaaaaa;   
                border-radius:3px;
                border-radius:3px;
                min-width:20;
            }
            QScrollBar::add-page:vertical,QScrollBar::sub-page:vertical,QScrollBar::add-page:horizontal,QScrollBar::sub-page:horizontal
            {
                background:#eeeeee;
                border-radius:4px;
            }
            """
# 导航条的样式
navi_str = """
            QWidget {
                border: none;
                background: #edeef0;
                color:#3265f6;
                padding:8px;
                font-size:8px} 
            QWidget:hover {
                background-color: #dddddd;
            }
            QWidget:pressed {
                background-color: #dddddd;
            }
            """
fig_navi_str = """
            QWidget {
                border: none;
                background: white;
                color:#3265f6;
                padding:2px;
                font-size:6px} 
            QWidget:hover {
                background-color: #dddddd;
            }
            QWidget:pressed {
                background-color: #dddddd;
            }
                
                """
# label的分子式的样式
label_struct_str = """QPushButton {
                        background-color:none;
                        color :black;
                        text-align:left;
                        padding: 1px 0px;
                        margin: 0px;
                    }
                    """
# label的分子式的样式被选中的样式
push_label_struct_str = """QPushButton {
                            background-color:none;
                            color :blue;
                            text-align:left;
                            padding: 1px 0px;
                            margin: 0px;
                        }
                        """

# 原始谱界面的 scan编号右侧的图标按钮
scan_icon_str = """
                QPushButton {
                    border-radius: 5px;
                    padding: 5px 5px;
                    background-color:none;
                    border-image:url(icon/QT1.png);
                    color: gray;
                }
                """

push_scan_icon_str = """
                QPushButton {
                    border-radius: 5px;
                    padding: 5px 5px;
                    background-color:none;
                    border-image:url(icon/QT.png);
                    color: gray;
                }
                QPushButton:hover {
                    border-image:url(icon/QT2.png);

                }
                QPushButton:pressed {
                     border-image:url(icon/QT2.png);
                }
                """
# 原始谱界面的 scan编号右侧的减号按钮被按下去时
scan_btn_str = """
                QPushButton {
                    border-radius: 4px;
                    padding: 5px 12px;
                    background-color: #eaeaea;
                    color: gray;
                }
                QPushButton:hover {
                    background-color: #dddddd;
                }
                QPushButton:pressed {
                    background-color: #dddddd;
                }
                """
