import pickle as pk
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import re

import time

MIN_SCAN_LEN = 5
MIN_MASS_LEN = 5
ZOOM_PARAM = 1.2
VIEW_LEN = 0.1
# [background, select, tic, axis_bg]
COLOR_LIST = ['skyblue', 'b', 'black']
COLOR_LIST1 = ['#B7B8B6', '#34675C', '#B3C100']
COLOR_LIST2 = ['#C1E1DC', '#265C00', '#68A225']
COLOR_LIST3 = ['#C1E1DC', '#217CA3', '#32384D']
COLOR_LIST4 = ['#B9C4C9', '#128277', '#004D47']
COLOR_LIST5 = ['#C1E1DC', '#2962FF', 'black', '#7b1e7a']
COLOR_LIST6 = ['#b9cfd4', '#7b1e7a', 'black', '#fffbff']
LINE_COLOR_LIST1 = ['#669900', '#99cc33', '#ccee66', '#006699', '#3399cc', '#990066', '#cc3399', '#ff6600', '#ff9900',
                    '#ffcc00']
LINE_COLOR_LIST2 = ['#ffbc42', '#d81159', '#8f2d56', '#218380', '#73d2de']
LINE_COLOR_LIST3 = ['#f6511d', '#ffb400', '#00a6ed', '#7fb800', '#0d2c54']
# LINE_COLOR_LIST4 = ['#fe5d26', '#f2c078', '#faedca', '#c1dbb3', '#7ebc89']
# LINE_COLOR_LIST_e = ['#6e44ff', '#b892ff', '#ffc2e2', '#ff90b3', '#ef7a85']
LINE_COLOR_LIST5 = ['#ffc107', '#007bff', '#17a2b8', '#28a745', '#dc3545']
LINE_COLOR_LIST6 = ['#e4572e', '#17bebb', '#ffc914', '#2e282a', '#76b041']


class EventFactory(object):
    def __init__(self, fig, flag):
        self.fig = fig  # fig 外层展示用
        self.active_flag = flag
        self.min_mass_len = MIN_MASS_LEN  # 放大时最小的mass跨度
        self.min_scan_len = MIN_SCAN_LEN  # 放大时最小的scan跨度
        self.zoom_param = ZOOM_PARAM  # 放大缩小灵敏参数
        self.view_len = VIEW_LEN  # 选中一个scan时左右view的长度参数
        self.line_color_list = COLOR_LIST6  # [background_line, select_line, tic_line, axis_bg]
        self.all_lines_color_list = LINE_COLOR_LIST5  # 所有线的颜色

    def init(self, ax, alldata, mass_range, scan_range):
        self.ax = ax  # 坐标轴
        self.ori_data = alldata  # 所有的数据

        self.colors_dict = {}  # 记录每条线的颜色，mass:color 避免出现放大缩小后同一条线颜色改变的问题
        self.ori_mass_range = mass_range
        self.ori_scan_range = scan_range  # 数据的mass和scan范围
        self.mass_range = self.ori_mass_range
        self.scan_range = self.ori_scan_range  # 画图的mass和scan范围
        self.draw_label = None  # 已经选中的类别
        self.draw_value = None  # 选中类别的值
        self.select_scan_range = []  # 选中的scan范围

        self.init_info()  # 初始化各类信息
        self.draw_all()  # 开始画图

    def init_info(self):
        self.last_position_xy = None  # 上次点击的位置
        self.last_press_time = None  # 上次点击的时间
        self.cuboid_lines = []  # 选择scan范围时，外框包络线*12， list of lines
        self.double_click_flag = False  # 双击的flag
        self.press_flag = False  # 鼠标按下的flag，用于选择scan范围
        self.scan_lines = []  # 选择scan时，着色的一段， list of lines
        self.line = None  # 选中的一条线
        self.label_move = None  # 移动时的label
        self.label_press = None  # 选择后的label
        self.points = []  # 选择scan的一个值时，标注连接点
        self.ori_len = None  # 当前图像上mass线的数目
        self.label_name = None  # 移动时标注，mass or scan
        self.value = None  # 与上面对应，具体的值是多少
        self.colors = []  # 当前mass线的颜色list， 用于着色
        self.tic_prop = None  # tic图缩放的比例

    def mouse_move(self, event):
        if self.active_flag[0] != 1:
            return
            # 鼠标移动事件
        if not event.inaxes:
            return
        self.del_text_while_move(event)
        label_name, value = self.resolve_location(event.xdata, event.ydata)
        if label_name is None:
            return
        if label_name == 'mass':
            tstr = 'mass = ' + str(value)
            self.label_move = event.inaxes.text3D(value,
                                                  -0.3 * (self.scan_range[1] - self.scan_range[0]) + self.scan_range[0],
                                                  0, tstr, bbox=dict(facecolor='green', alpha=0.3))
            self.label_name = 'mass'
            self.value = value
        else:
            tstr = 'scan = ' + str(value)
            self.label_move = event.inaxes.text3D((self.mass_range[1] - self.mass_range[0]) * 0.2 + self.mass_range[1],
                                                  value, 0, tstr, bbox=dict(facecolor='green', alpha=0.3))
            self.label_name = 'scan'
            self.value = value

        if self.press_flag and label_name == 'scan':
            self.draw_cuboid(value)

        self.set_lim()
        self.fig.canvas.draw_idle()

    def draw_cuboid(self, value):
        # 画正方体
        label1, value1 = self.resolve_location(self.last_position_xy[0], self.last_position_xy[1])
        if label1 != 'scan':
            return

        if abs(value - value1) < 1:
            return
        z_min = int(self.ax.get_zlim3d()[0])
        z_max = int(self.ax.get_zlim3d()[1])
        x_min = self.mass_range[0]
        x_max = self.mass_range[1]
        y1 = value
        y2 = value1

        self.cuboid_lines[0].set_data_3d(np.array([x_min, x_min]), np.array([y1, y2]), np.array([z_min, z_min]))
        self.cuboid_lines[1].set_data_3d(np.array([x_min, x_max]), np.array([y2, y2]), np.array([z_min, z_min]))
        self.cuboid_lines[2].set_data_3d(np.array([x_max, x_max]), np.array([y2, y1]), np.array([z_min, z_min]))
        self.cuboid_lines[3].set_data_3d(np.array([x_max, x_min]), np.array([y1, y1]), np.array([z_min, z_min]))
        self.cuboid_lines[4].set_data_3d(np.array([x_min, x_min]), np.array([y1, y1]), np.array([z_min, z_max]))
        self.cuboid_lines[5].set_data_3d(np.array([x_min, x_min]), np.array([y2, y2]), np.array([z_min, z_max]))
        self.cuboid_lines[6].set_data_3d(np.array([x_max, x_max]), np.array([y1, y1]), np.array([z_min, z_max]))
        self.cuboid_lines[7].set_data_3d(np.array([x_max, x_max]), np.array([y2, y2]), np.array([z_min, z_max]))
        self.cuboid_lines[8].set_data_3d(np.array([x_min, x_min]), np.array([y1, y2]), np.array([z_max, z_max]))
        self.cuboid_lines[9].set_data_3d(np.array([x_min, x_max]), np.array([y2, y2]), np.array([z_max, z_max]))
        self.cuboid_lines[10].set_data_3d(np.array([x_max, x_max]), np.array([y2, y1]), np.array([z_max, z_max]))
        self.cuboid_lines[11].set_data_3d(np.array([x_max, x_min]), np.array([y1, y1]), np.array([z_max, z_max]))

    def mouse_press(self, event):
        if self.active_flag[0] != 1:
            return
            # 鼠标点击事件
        self.press_flag = False
        if not event.inaxes:
            self.last_position_xy = None
            self.last_press_time = None
            return

        this_press_time = time.time()
        if self.last_press_time is None:
            self.last_press_time = this_press_time
            self.last_position_xy = (event.xdata, event.ydata)
            self.press_flag = True

        elif (this_press_time - self.last_press_time) < 0.5 and self.dis(self.last_position_xy,
                                                                         (event.xdata, event.ydata)) < 10E-3:
            self.double_click_flag = True
            self.double_click()
            return

        else:
            self.last_press_time = this_press_time
            self.last_position_xy = (event.xdata, event.ydata)
            self.press_flag = True

    def mouse_release(self, event):
        if self.active_flag[0] != 1:
            return
        # 鼠标释放事件
        self.press_flag = False
        if not event.inaxes:
            self.last_position_xy = None
            self.last_press_time = None
            return
        release_time = time.time()
        if self.double_click_flag:
            self.double_click_flag = False
            return
        if self.last_press_time is None:
            return
        if release_time - self.last_press_time < 0.3 and self.dis(self.last_position_xy,
                                                                  (event.xdata, event.ydata)) < 10E-5:
            self.click_on(event)
        elif release_time - self.last_press_time > 0.5 and self.dis(self.last_position_xy,
                                                                    (event.xdata, event.ydata)) > 10E-5:
            # self.zoom_by_mouse(event)
            self.select_by_mouse(event)

    def scroll(self, event):
        if self.active_flag[0] != 1:
            return
        # 滚轮事件
        if not event.inaxes:
            return
        self.zoom_by_scroll(event)

    def click_on(self, event):
        # 单击事件
        self.del_line_while_click()
        if self.label_name is None:
            return
        elif self.label_name == 'mass':
            self.draw_mass()

        else:
            self.draw_scan()

        self.set_lim()
        self.fig.canvas.draw_idle()

    def double_click(self):
        # 双击事件
        self.draw_label = None
        self.draw_value = None
        self.select_scan_range = []
        self.mass_range = self.ori_mass_range
        self.scan_range = self.ori_scan_range

        self.draw_all()
        self.double_click_flag = True

    def draw_all(self):
        ## 画整张图，初始化其他状态
        t = time.time()
        self.all_data = [[], [], []]
        # 准备数据
        for i, mass in enumerate(self.ori_data[0]):
            if mass < self.mass_range[0] or mass > self.mass_range[1]:
                continue
            self.all_data[0].append(mass)
            self.all_data[2].append(self.ori_data[2][i][self.scan_range[0]:self.scan_range[1]])

        self.all_data[1] = self.ori_data[1][self.scan_range[0]:self.scan_range[1]]

        # 清空画板,初始化各个状态
        self.ax.cla()
        self.init_info()
        self.ax.xaxis.pane.fill = True
        self.ax.yaxis.pane.fill = True
        self.ax.zaxis.pane.fill = True
        self.ax.xaxis.pane.set_color(self.line_color_list[3])
        self.ax.yaxis.pane.set_color(self.line_color_list[3])
        self.ax.zaxis.pane.set_color(self.line_color_list[3])
        self.mass_lim = 0.1 * (self.mass_range[1] - self.mass_range[0])
        self.scan_lim = 0.1 * (self.scan_range[1] - self.scan_range[0])
        view_len = int(self.view_len * (self.scan_range[1] - self.scan_range[0]))
        self.view_range = [-view_len, view_len]

        # 画上原始线

        self.ax.set_xlabel(' mass ', labelpad=15)
        self.ax.set_ylabel(' scan ', labelpad=15)
        self.ax.set_zlabel('intensity', labelpad=20)
        x = self.all_data[0]
        y = self.all_data[1]
        z = self.all_data[2]

        tmp_y = np.array(y)
        ylimit = len(y)
        L = len(x)
        max_z_list = []
        len_color = len(self.all_lines_color_list)
        for i in range(L):
            tmp_x = np.full(ylimit, x[i])
            tmp_z = np.array(z[i])
            max_z_list.append(max(z[i]))
            self.ax.plot(tmp_x, tmp_y, tmp_z, linewidth=0.4, zorder=1, color=self.all_lines_color_list[i % len_color])
            # self.ax.plot(tmp_x, tmp_y, tmp_z, linewidth=0.4, zorder=1)

        max_z = max(max_z_list)
        if len(self.colors_dict) == 0:
            for line in self.ax.lines:
                self.colors_dict[line.get_data_3d()[0][0]] = line.get_color()
        else:
            for line in self.ax.lines:
                line.set_color(self.colors_dict[line.get_data_3d()[0][0]])

        self.ori_len = len(self.ax.lines)
        for i in range(self.ori_len):
            tline = \
            self.ax.plot(np.array([0]), np.array([0]), np.array([0]), color=self.ax.lines[i].get_color(), linewidth=0.6,
                         zorder=2)[0]
            self.scan_lines.append(tline)

        self.line = \
        self.ax.plot(np.array([0]), np.array([0]), np.array([0]), color=self.line_color_list[1], linewidth=1.0,
                     zorder=3)[0]
        # 画出TIC
        max_tic = max(self.ori_data[3][self.scan_range[0]:self.scan_range[1]])
        self.tic_prop = max_z / max_tic
        tmp_z = np.array(self.ori_data[3][self.scan_range[0]:self.scan_range[1]]) * self.tic_prop
        tmp_x = np.full(ylimit, self.mass_range[0])
        self.tic_line = self.ax.plot(tmp_x, tmp_y, tmp_z, color=self.line_color_list[2], linewidth=0.6, zorder=0)[0]
        self.scan_tic = \
        self.ax.plot(np.array([0]), np.array([0]), np.array([0]), linestyle="-.", color=self.line_color_list[2],
                     linewidth=0.4, zorder=0)[0]

        self.cuboid_lines = []
        for i in range(12):
            tmp_line = \
            self.ax.plot(np.array([0]), np.array([0]), np.array([0]), linestyle="-.", color='black', linewidth=0.8,
                         zorder=4)[0]
            self.cuboid_lines.append(tmp_line)
        for line in self.ax.lines:
            self.colors.append(line.get_color())
        # 画一个长方体的12条棱，用于拖拽选择多个scan时的对齐

        if self.draw_label == 'mass' and self.mass_range[0] < self.draw_value < self.mass_range[1]:
            self.value = self.draw_value
            self.label_name = self.draw_label
            self.draw_mass()
        elif self.draw_label == 'scan' and self.scan_range[0] < self.draw_value < self.scan_range[1]:
            self.value = self.draw_value
            self.label_name = self.draw_label
            self.draw_scan()

        elif self.draw_label == 'scan_range' and (
                self.scan_range[0] < self.select_scan_range[0] < self.scan_range[1] or self.scan_range[0] <
                self.select_scan_range[1] < self.scan_range[1]):
            self.draw_scan_range()

        elif self.draw_label == 'mass' or self.draw_label == 'scan' or self.draw_label == 'scan_range':
            for i in range(self.ori_len):
                self.ax.lines[i].set_color(self.line_color_list[0])

        self.ax.zaxis.get_major_formatter().set_powerlimits((0, 1))
        self.ax.zaxis.get_major_formatter().set_useMathText(True)
        self.ax.tick_params(axis='y', labelrotation=30)
        self.ax.tick_params(axis='x', labelrotation=-10)
        self.ax.tick_params(axis='z', labelrotation=-3.5)
        self.set_lim()
        self.fig.canvas.draw_idle()
        print(time.time() - t)

    def draw_mass(self):
        # 在mass上选择一条线，画出该线
        y_data = None
        z_data = None
        for line in self.ax.lines:
            tmp = line.get_data_3d()[0][0]
            if tmp == self.value:
                y_data = line.get_data_3d()[1]
                z_data = line.get_data_3d()[2]
                x_data = np.full(len(z_data), self.value)
                break

        self.line.set_data_3d(x_data, y_data, z_data)
        tstr = self.label_name + ' = ' + str(self.value)
        self.label_press = self.ax.text3D(self.value,
                                          -0.3 * (self.scan_range[1] - self.scan_range[0]) + self.scan_range[0], 0,
                                          tstr, bbox=dict(facecolor='red', alpha=0.3))
        for i in range(self.ori_len):
            self.ax.lines[i].set_color(self.line_color_list[0])

        self.line.set_color(self.line_color_list[1])
        self.draw_label = 'mass'
        self.draw_value = self.value

    def draw_scan(self):
        # 在scan上选择一根线，画出该线
        x_data = []
        y_data = []
        z_data = []
        begin = self.value + self.view_range[0]
        if begin < self.scan_range[0]:
            begin = self.scan_range[0]

        end = self.value + self.view_range[1]
        if end > self.scan_range[1]:
            end = self.scan_range[1]

        L = end - begin
        x_data.append(self.mass_range[0])
        y_data.append(self.value)
        z_data.append(0)
        for i, item in enumerate(self.all_data[0]):
            x_data.append(item)
            x_data.append(item)
            x_data.append(item)
            y_data.append(self.value)
            y_data.append(self.value)
            y_data.append(self.value)
            z_data.append(0)
            z_data.append(self.all_data[2][i][self.value - self.scan_range[0]])

            z_data.append(0)
            tmp_x = np.full(L, item)
            tmp_y = np.array(list(range(begin, end)))
            tmp_z = np.array(self.all_data[2][i][begin - self.scan_range[0]:end - self.scan_range[0]])

            self.scan_lines[i].set_data_3d(tmp_x, tmp_y, tmp_z)
            self.points.append(self.ax.scatter(item, self.value, self.all_data[2][i][self.value - self.scan_range[0]],
                                               c=self.colors[i], s=8, marker='_'))

        x_data.append(self.mass_range[1])
        y_data.append(self.value)
        z_data.append(0)
        x_data = np.array(x_data)
        y_data = np.array(y_data)
        z_data = np.array(z_data)
        self.line.set_data_3d(x_data, y_data, z_data)
        self.scan_tic.set_data_3d(np.array([self.mass_range[0], self.mass_range[0]]),
                                  np.array([self.value, self.value]),
                                  np.array([0, self.ori_data[3][self.value] * self.tic_prop]))
        tstr = self.label_name + ' = ' + str(self.value)
        self.label_press = self.ax.text3D((self.mass_range[1] - self.mass_range[0]) * 0.2 + self.mass_range[1],
                                          self.value, 0, tstr, bbox=dict(facecolor='red', alpha=0.3))
        for i in range(self.ori_len):
            self.ax.lines[i].set_color(self.line_color_list[0])

        self.draw_label = 'scan'
        self.draw_value = self.value

    def draw_scan_range(self):
        # 在scan上选择一段，画出该段
        begin = max(self.select_scan_range[0], self.scan_range[0])
        end = min(self.select_scan_range[1], self.scan_range[1])

        L = end - begin
        for i, item in enumerate(self.all_data[0]):
            tmp_x = np.full(L, item)
            tmp_y = np.array(list(range(begin, end)))
            tmp_z = np.array(self.all_data[2][i][begin - self.scan_range[0]:end - self.scan_range[0]])
            self.scan_lines[i].set_data_3d(tmp_x, tmp_y, tmp_z)

        self.scan_tic.set_data_3d(
            np.array([self.mass_range[0], self.mass_range[0], self.mass_range[0], self.mass_range[0]]),
            np.array([begin, begin, end, end]),
            np.array([self.ori_data[3][begin] * self.tic_prop, 0, 0, self.ori_data[3][end] * self.tic_prop]))
        tstr = 'scan:' + str(begin) + '-' + str(end)
        self.label_press = self.ax.text3D((self.mass_range[1] - self.mass_range[0]) * 0.2 + self.mass_range[1],
                                          int((begin + end) / 2), 0, tstr, bbox=dict(facecolor='red', alpha=0.3))
        for i in range(self.ori_len):
            self.ax.lines[i].set_color(self.line_color_list[0])

        self.draw_label = 'scan_range'
        self.draw_value = None

    def zoom_by_mouse(self, event):
        # 利用鼠标拖拽 放大或者缩小，当前不可用
        if self.last_position_xy is None:
            return
        label1, value1 = self.resolve_location(self.last_position_xy[0], self.last_position_xy[1])
        if label1 is None:
            return
        label2, value2 = self.resolve_location(event.xdata, event.ydata)
        if label2 != label1:
            return

        tstr1 = self.ax.format_coord(self.last_position_xy[0], self.last_position_xy[1])
        tstr2 = self.ax.format_coord(event.xdata, event.ydata)

        ans1 = re.match('x=(.*) , y=(.*), z=(.*)', tstr1)
        ans2 = re.match('x=(.*) , y=(.*), z=(.*)', tstr2)
        if ans1 is None or ans2 is None:
            return

        if label1 == 'mass':
            x1 = float(ans1.group(1))
            x2 = float(ans2.group(1))
            x1 = int(x1)
            x2 = int(x2)

            self.mass_range = [min(x1, x2), max(x1, x2)]

        elif label1 == 'scan':
            y1 = float(ans1.group(2))
            y2 = float(ans2.group(2))
            y1 = int(y1)
            y2 = int(y2)

            self.scan_range = [min(y1, y2), max(y1, y2)]

        self.draw_all()

    def select_by_mouse(self, event):
        # 利用鼠标拖拽，选择一段
        self.del_line_while_click()
        if self.last_position_xy is None:
            return
        label1, value1 = self.resolve_location(self.last_position_xy[0], self.last_position_xy[1])
        if label1 is None:
            return
        label2, value2 = self.resolve_location(event.xdata, event.ydata)
        if label2 != label1:
            return

        if label2 != 'scan':
            return

        tstr1 = self.ax.format_coord(self.last_position_xy[0], self.last_position_xy[1])
        tstr2 = self.ax.format_coord(event.xdata, event.ydata)

        ans1 = re.match('x=(.*) , y=(.*), z=(.*)', tstr1)
        ans2 = re.match('x=(.*) , y=(.*), z=(.*)', tstr2)
        if ans1 is None or ans2 is None:
            return

        y1 = float(ans1.group(2))
        y2 = float(ans2.group(2))
        y1 = int(y1)
        y2 = int(y2)

        self.select_scan_range = [min(y1, y2), max(y1, y2)]
        # 在图上表现出来
        self.draw_scan_range()
        self.last_position_xy = None
        self.last_press_time = None
        self.double_click_flag = False
        self.set_lim()
        self.fig.canvas.draw_idle()

    def zoom_by_scroll(self, event):
        # 利用滚轮完成放大缩小
        label, value = self.resolve_location(event.xdata, event.ydata)
        if label is None:
            return
        tstr = self.ax.format_coord(event.xdata, event.ydata)
        ans = re.match('x=(.*) , y=(.*), z=(.*)', tstr)

        if label == 'mass':
            low, high = self.calculation_range(self.mass_range, int(float(ans.group(1))), event.button,
                                               self.min_mass_len, self.ori_mass_range)
            self.mass_range = [low, high]
            if self.mass_range == self.ori_mass_range:
                return

        elif label == 'scan':
            low, high = self.calculation_range(self.scan_range, int(float(ans.group(2))), event.button,
                                               self.min_scan_len, self.ori_scan_range)
            self.scan_range = [low, high]
            if self.scan_range == self.ori_scan_range:
                return

        self.draw_all()

    def del_text_while_move(self, event):
        # 移动鼠标时，实时更新label信息时，删除原有label信息
        if self.label_move:
            self.label_move.set_text('')
            event.inaxes.texts.remove(self.label_move)

        for line in self.cuboid_lines:
            line.set_data_3d(np.array([0]), np.array([0]), np.array([0]))

        self.label_move = None
        self.label_name = None
        self.value = None
        self.set_lim()
        self.fig.canvas.draw_idle()

    def del_line_while_click(self):
        # 鼠标单击or拖拽时，删除原有信息
        if self.label_press:
            self.label_press.set_text('')
            self.ax.texts.remove(self.label_press)
        L = len(self.points)
        for i in range(L):
            self.points[i].remove()

        for i in range(L):
            item = self.points.pop()
            del (item)

        for line in self.scan_lines:
            line.set_data_3d(np.array([0]), np.array([0]), np.array([0]))
        self.scan_tic.set_data_3d(np.array([0]), np.array([0]), np.array([0]))
        self.label_press = None
        self.draw_label = None
        self.draw_value = None
        self.select_scan_range = []
        i = 0
        for item in self.ax.lines:
            item.set_color(self.colors[i])
            i += 1

        self.line.set_data_3d(np.array([0]), np.array([0]), np.array([0]))
        self.set_lim()
        self.fig.canvas.draw_idle()

    def resolve_location(self, xdata, ydata):
        # 判断落点的位置是否符合要求，返回落点属于mass or scan以及对应的值
        tstr = self.ax.format_coord(xdata, ydata)
        ans = re.match('x=(.*) , y=(.*), z=(.*)', tstr)
        if ans is None:
            return None, None

        x = float(ans.group(1))
        y = float(ans.group(2))
        z = float(ans.group(3))

        A = ((self.mass_range[0] - 0.5) <= x <= self.mass_range[1]) and (
                self.scan_range[0] - self.scan_lim <= y <= self.scan_range[0] + self.scan_lim)
        B = ((self.mass_range[1] - self.mass_lim) <= x <= (self.mass_range[1] + self.mass_lim)) and (
                self.scan_range[0] <= y <= (self.scan_range[1] + 0.5))

        if not A and not B:
            return None, None

        label_name = None
        if B:
            value = int(round(y))
            if value < self.scan_range[0]:
                value = self.scan_range[0]
            if value > self.scan_range[1]:
                value = self.scan_range[1]
            return 'scan', value

        else:
            label_name = 'mass'
            if x > self.mass_range[1] + 1 or x < self.mass_range[0] - 1:
                return None, None
            minm = self.mass_range[1]
            mass = 0
            for line in self.ax.lines:
                tmp = line.get_data_3d()[0][0]
                if abs(tmp - x) < minm:
                    minm = abs(tmp - x)
                    mass = tmp

            return 'mass', mass

    def calculation_range(self, ori_range, aim_pos, button, min_len, max_range):
        # 计算缩放的范围
        ori_len = ori_range[1] - ori_range[0] + 1
        if button == 'down':
            # 放大
            aim_len = int(ori_len * self.zoom_param)
        else:
            # 缩小
            aim_len = int(ori_len / self.zoom_param)

        if aim_len >= (max_range[1] - max_range[0] + 1):
            return max_range[0], max_range[1]

        elif aim_len <= min_len:
            aim_len = min_len

        left = int(((aim_pos - ori_range[0] + 1) / ori_len) * aim_len)
        right = aim_len - left
        low = aim_pos - left
        if low < max_range[0]:
            tmp = max_range[0] - low
            low = max_range[0]
            right = right + tmp

        high = aim_pos + right
        if high > max_range[1]:
            tmp = high - max_range[1]
            high = max_range[1]
            left = max(left - tmp, max_range[0])

        return low, high

    def dis(self, pos1, pos2):
        # 计算鼠标两次点击位置的距离
        x = (pos1[0] - pos2[0]) ** 2
        y = (pos1[1] - pos2[1]) ** 2
        ret = np.sqrt(x + y)
        return ret

    def set_lim(self):
        # 设置坐标轴的limit
        self.ax.set_xlim3d(left=self.mass_range[0])
        self.ax.set_xlim3d(right=self.mass_range[1])
        self.ax.set_ylim3d(bottom=self.scan_range[0])
        self.ax.set_ylim3d(top=self.scan_range[1])
        self.ax.set_zlim3d(bottom=0)

    def get_current_scan(self):
        # 外部事件，获取当前选择的scan值
        if self.draw_label == 'scan':
            return self.draw_value
        # 考虑选择的线不在当前范围内
        return None

    def get_current_scan_range(self):
        # 外部事件，获取当前选择的scan范围
        if len(self.select_scan_range) == 2:
            return self.select_scan_range
        return None
        # 考虑选择的scan不在当前选择的范围内
