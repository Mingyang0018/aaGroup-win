"""
@File : aaGroup.py
@Time : 2023/5/20 23:03:01
@Author : <Yang Mingxuan>
@Version : 1.0
@description: 氨基酸序列同源聚类分析
@update: 
"""

# 导入必要的模块
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from Bio import Align
from multiprocessing import Pool
import mmap,os,dill,sys,shutil,pathlib
from pathlib import Path
from PyQt5.QtWidgets import QApplication, QWidget, QFileDialog, QLabel, QLineEdit, QPushButton, QComboBox, QMessageBox, QMainWindow, QAction, QTableWidget,QTableWidgetItem,QSizePolicy,QGridLayout,QGroupBox,QHBoxLayout,QTextBrowser,QFrame,QTabBar,QStackedWidget
from PyQt5.QtGui import QPixmap, QDesktopServices
from PyQt5.QtCore import Qt,QUrl, pyqtSignal
import datetime

from PyQt5.QtWebEngineWidgets import QWebEngineView
import ydata_profiling

import tempfile
import seaborn as sns

import matplotlib,random 
# matplotlib.use("Qt5Agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg 
from sklearn.cluster import KMeans,MiniBatchKMeans
import multiprocessing

EXCEL_TYPE = ["xls", "xlsx", "xlsm","csv"]

# 定义一个继承 QLabel 的类
class ClickableLabel(QLabel):
    # 定义一个自定义的 clicked 信号，可以传递参数
    clicked = pyqtSignal(str)

    # 重写 mouseReleaseEvent 函数
    def mouseReleaseEvent(self, event):
        # 调用父类的函数
        super().mouseReleaseEvent(event)
        # 如果是左键释放，就发出 clicked 信号，并传递标签的文本
        if event.button() == Qt.LeftButton:
            self.clicked.emit(self.text())

# 定义一个前端类
class MainView(QMainWindow):
    def __init__(self):
        super().__init__()
        self.abs_path = os.path.abspath(__file__)
        self.parent_path = os.path.dirname(self.abs_path)
        self.initUI()

    def initUI(self):
        # 设置窗口标题和大小
        self.setWindowTitle('氨基酸同源分析器1.0')
        self.resize(410,350)
        
        # 创建一个菜单栏对象
        self.menuBar = self.menuBar()
        # 添加一个帮助菜单，并返回一个 QMenu 对象
        self.helpMenu = self.menuBar.addMenu('帮助')
        # 创建一个 QAction 对象，并设置文本和快捷键
        self.helpAction = QAction('查看帮助', self)
        # 把 QAction 对象添加到 QMenu 对象中
        self.helpMenu.addAction(self.helpAction)
        # 把 QAction 对象的 triggered 信号，连接到 show_help 方法上
        self.helpAction.triggered.connect(self.show_help)
        # 创建中央部件和布局
        self.centralWidget = QWidget(self)
        self.setCentralWidget(self.centralWidget)
        self.layout = QGridLayout(self.centralWidget)
        self.layout.setSpacing(10)
        # 创建标签栏和堆叠部件
        self.tabBar = QTabBar(self.centralWidget)
        # 添加标签项
        self.tabBar.addTab("文件合并")
        self.tabBar.addTab("同源分组")
        self.tabBar.addTab("数据分析")
    
        # 设置标签项的样式表
        self.tabBar.setStyleSheet("""
            QTabBar::tab {
                background: lightgray;
                border: none;
                padding: 7px;
            }
            QTabBar::tab:selected {
                background: white;
            }
            """)
        # 设置标签栏的尺寸策略
        self.tabBar.setExpanding(True) # 横向充满整个窗口宽度
        # 将标签栏添加到布局中
        self.layout.addWidget(self.tabBar, 0, 0, 1, 3) # 占据三列
        # 创建堆叠部件和子窗口部件
        self.stackedWidget = QStackedWidget(self.centralWidget)
        # 调用函数创建子窗口部件和标签
        self.createSubWindow(model=0)
        self.createSubWindow(model=1)
        self.createSubWindow(model=2)
        # 将堆叠部件添加到布局中
        self.layout.addWidget(self.stackedWidget, 1, 0, 1, 3) # 占据三列
        # 设置初始显示的子窗口
        self.stackedWidget.setCurrentIndex(0)
        # 连接信号和槽函数
        self.tabBar.currentChanged.connect(self.switchToSubWidget)

    # 定义显示子窗口的槽函数
    def switchToSubWidget(self, index):
        # 根据索引切换子窗口
        self.stackedWidget.setCurrentIndex(index)

    # 定义创建子窗口部件和标签的函数
    def createSubWindow(self, model):
        self.model = model
        layout = self.layout
        if self.model == 0:
            # 创建子窗口部件和标签
            self.subWindow01 = QWidget(self.stackedWidget)
            self.subLayout01 = QGridLayout(self.subWindow01)
            self.subLayout01.setColumnStretch(0,10)
            self.subLayout01.setColumnStretch(1,10)
            self.subLayout01.setColumnStretch(2,10)
            # 创建标签和输入框
            self.label01_1 = QLabel('原始文件夹:', self.subWindow01)
            # self.label01_1.setAlignment(Qt.AlignCenter)
            self.edit01_1 = QLineEdit(self.subWindow01)
            
            # 创建按钮和绑定事件
            self.button01_1 = QPushButton('选择路径', self.subWindow01)
            self.button01_1.clicked.connect(self.select_source_folder)
            
            self.label01_2 = QLabel('是否去重:', self.subWindow01)
            self.combo01_2 = QComboBox(self)
            self.combo01_2.addItems(["是", '否'])

            self.label01_3 = QLabel('去重列名:', self.subWindow01)
            self.edit01_3 = QLineEdit(self)
            self.edit01_4 = QLineEdit(self)

            # self.label4 = QLabel('感谢老板打赏!', self.subWindow02)
            self.labelimage01_4 = QLabel('', self.subWindow01)
            pixmap = QPixmap(self.parent_path+"/GZlogo.jpg")
            self.labelimage01_4.setPixmap(pixmap.scaled(120, 120))

            # 创建下载地址链接
            self.url_label01_4 = ClickableLabel(self.subWindow01)
            self.url_label01_4.setText('<font color="blue">下载地址</font>')
            # 把网址作为一个属性存储在 self.url_label 中
            self.url_label01_4.url = 'https://pan.baidu.com/s/1BckIYp4DQvonMBAHFFzrJw?pwd=1234'
            # 设置鼠标样式为手型
            self.url_label01_4.setCursor(Qt.PointingHandCursor)
            # 连接 clicked 信号和 open_url 槽函数
            self.url_label01_4.clicked.connect(lambda: self.open_url(self.url_label01_4.url))
            
            self.button01_2 = QPushButton('文件合并', self.subWindow01)
            self.button01_2.clicked.connect(lambda:self.click(0))

            self.subLayout01.addWidget(self.label01_1,1,0,1,1) # 将标签添加到网格布局中，占用第二行第一列
            self.subLayout01.addWidget(self.edit01_1,1,1,1,1) # 将输入框添加到网格布局中，占用第二行第二列
            self.subLayout01.addWidget(self.button01_1,1,2,1,1)
            self.subLayout01.addWidget(self.label01_2,2,0,1,1)
            self.subLayout01.addWidget(self.combo01_2,2,1,1,1)
            self.subLayout01.addWidget(self.label01_3,3,0,1,1)
            self.subLayout01.addWidget(self.edit01_3,3,1,1,1)
            self.subLayout01.addWidget(self.edit01_4,3,2,1,1)
            self.subLayout01.addWidget(self.labelimage01_4,4,1,2,2)
            self.subLayout01.addWidget(self.url_label01_4,4,0,1,1)
            self.subLayout01.addWidget(self.button01_2,7,1,1,1)

            self.stackedWidget.addWidget(self.subWindow01)

            # 创建子窗口对象
            self.subWindow01_sub_window = QWidget()
            self.subWindow01_sub_window.setWindowTitle("文件预览")
            self.subWindow01_sub_window.resize(750, 500)

            # 创建子窗口中的表格对象
            self.table = QTableWidget(self.subWindow01_sub_window)
            self.table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            # 设置表格的滚动条策略
            self.table.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
            self.table.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
            # 创建子窗口中的布局对象
            self.layout_subWindow01_sub_window = QGridLayout(self.subWindow01_sub_window)
            # 设置布局中的边距为0
            self.layout_subWindow01_sub_window.setContentsMargins(0, 0, 0, 50)
            
            # 创建子窗口中的确认按钮对象
            self.confirm_button01 = QPushButton("导出", self.subWindow01_sub_window)
            self.confirm_button01.clicked.connect(self.export_excel)
            self.layout_subWindow01_sub_window.addWidget(self.confirm_button01,0,0,1,1)
            # 将表格添加到布局中
            self.layout_subWindow01_sub_window.addWidget(self.table,1,0,1,10)
            
        elif self.model == 1:
            # 创建子窗口部件和标签
            self.subWindow02 = QWidget(self.stackedWidget)
            self.subLayout02 = QGridLayout(self.subWindow02)
            self.subLayout02.setColumnStretch(0,10)
            self.subLayout02.setColumnStretch(1,10)
            self.subLayout02.setColumnStretch(2,10)

            # 设置第二组界面
            # 创建标签和输入框
            self.label02_1 = QLabel('数据文件:', self.subWindow02)
            self.edit02_1 = QLineEdit(self.subWindow02)
            # 创建按钮和绑定事件
            self.button02_1 = QPushButton('选择文件', self.subWindow02)
            self.button02_1.clicked.connect(self.select_file)

            self.label02_2 = QLabel('类别:', self.subWindow02)
            self.combo02_2 = QComboBox(self.subWindow02)
            self.combo02_2.addItems(["无"])

            self.label02_3 = QLabel('氨基酸序列:', self.subWindow02)
            self.combo02_3 = QComboBox(self.subWindow02)
            # self.combo02_3.addItems(["无"])
            
            self.label02_4 = QLabel('分组数量:', self.subWindow02)
            self.edit02_4 = QLineEdit(self.subWindow02)
            self.edit02_4.setPlaceholderText("输入整数")
            
            # self.label4 = QLabel('感谢老板打赏!', self.subWindow02)
            self.labelimage02_5 = QLabel('', self.subWindow02)
            pixmap = QPixmap(self.parent_path+"/GZlogo.jpg")
            self.labelimage02_5.setPixmap(pixmap.scaled(120, 120))

            # 创建下载地址链接
            self.url_label02_5 = ClickableLabel(self.subWindow02)
            self.url_label02_5.setText('<font color="blue">下载地址</font>')
            # 把网址作为一个属性存储在 self.url_label 中
            self.url_label02_5.url = 'https://pan.baidu.com/s/1BckIYp4DQvonMBAHFFzrJw?pwd=1234'
            # 设置鼠标样式为手型
            self.url_label02_5.setCursor(Qt.PointingHandCursor)
            # 连接 clicked 信号和 open_url 槽函数
            self.url_label02_5.clicked.connect(lambda: self.open_url(self.url_label02_5.url))
            
            self.button02_6 = QPushButton('同源分组', self.subWindow02)
            self.button02_6.clicked.connect(lambda:self.click(1))

            self.subLayout02.addWidget(self.label02_1,1,0,1,1)
            self.subLayout02.addWidget(self.edit02_1,1,1,1,1)
            self.subLayout02.addWidget(self.button02_1,1,2,1,1)
            self.subLayout02.addWidget(self.label02_2,2,0,1,1)
            self.subLayout02.addWidget(self.combo02_2,2,1,1,1)
            self.subLayout02.addWidget(self.label02_3,3,0,1,1)
            self.subLayout02.addWidget(self.combo02_3,3,1,1,1)
            self.subLayout02.addWidget(self.label02_4,4,0,1,1)
            self.subLayout02.addWidget(self.edit02_4,4,1,1,1)
            # self.subLayout02.addWidget(self.label4,5,0,1,1)
            self.subLayout02.addWidget(self.labelimage02_5,5,1,2,2)
            self.subLayout02.addWidget(self.url_label02_5,5,0,1,1)
            self.subLayout02.addWidget(self.button02_6,8,1,1,1)

            self.stackedWidget.addWidget(self.subWindow02)
        elif self.model == 2:
            # 创建子窗口部件和标签
            self.subWindow03 = QWidget(self.stackedWidget)
            self.subLayout03 = QGridLayout(self.subWindow03)
            self.subLayout03.setColumnStretch(0,10)
            self.subLayout03.setColumnStretch(1,10)
            self.subLayout03.setColumnStretch(2,10)

            # 设置第二组界面
            # 创建标签和输入框
            self.label03_1 = QLabel('数据文件:', self.subWindow03)
            self.edit03_1 = QLineEdit(self.subWindow03)
            # 创建按钮和绑定事件
            self.button03_1 = QPushButton('选择文件', self.subWindow03)
            self.button03_1.clicked.connect(self.select_file02)

            self.label03_2 = QLabel('缺失值处理:', self)
            self.combo03_2 = QComboBox(self)
            self.combo03_2.addItems(
                ["无", "去除缺失值", '填充固定值', '填充均值', '填充中位数', '填充前值', '填充后值'])
            self.edit03_2 = QLineEdit(self)
            self.edit03_2.setPlaceholderText("输入固定值")
 
            self.label03_3 = QLabel('时间序列:', self)
            self.combo03_3 = QComboBox(self)
            self.combo03_3.addItems(["无", "智能识别"])

            self.labelimage03_5 = QLabel('', self.subWindow03)
            pixmap = QPixmap(self.parent_path+"/GZlogo.jpg")
            self.labelimage03_5.setPixmap(pixmap.scaled(120, 120))

            # 创建下载地址链接
            self.url_label03_5 = ClickableLabel(self.subWindow03)
            self.url_label03_5.setText('<font color="blue">下载地址</font>')
            # 把网址作为一个属性存储在 self.url_label 中
            self.url_label03_5.url = 'https://pan.baidu.com/s/1BckIYp4DQvonMBAHFFzrJw?pwd=1234'
            # 设置鼠标样式为手型
            self.url_label03_5.setCursor(Qt.PointingHandCursor)
            # 连接 clicked 信号和 open_url 槽函数
            self.url_label03_5.clicked.connect(lambda: self.open_url(self.url_label03_5.url))
            
            self.button03_6 = QPushButton('数据分析', self.subWindow03)
            self.button03_6.clicked.connect(lambda:self.click(2))

            self.subLayout03.addWidget(self.label03_1,1,0,1,1)
            self.subLayout03.addWidget(self.edit03_1,1,1,1,1)
            self.subLayout03.addWidget(self.button03_1,1,2,1,1)
            self.subLayout03.addWidget(self.label03_2,2,0,1,1)
            self.subLayout03.addWidget(self.combo03_2,2,1,1,1)
            self.subLayout03.addWidget(self.edit03_2,2,2,1,1)
            self.subLayout03.addWidget(self.label03_3,3,0,1,1)
            self.subLayout03.addWidget(self.combo03_3,3,1,1,1)
            self.subLayout03.addWidget(self.labelimage03_5,5,1,2,2)
            self.subLayout03.addWidget(self.url_label03_5,5,0,1,1)
            self.subLayout03.addWidget(self.button03_6,8,1,1,1)

            self.stackedWidget.addWidget(self.subWindow03)
        else:
            pass

    def select_file02(self):
        # 弹出一个文件对话框，返回文件路径
        file_path, _ = QFileDialog.getOpenFileName(self, '选择文件')
        if file_path:
            # self.file_path = file_path
            self.edit03_1.setText(file_path)
            file_ext = file_path.split('.')[-1]  # 获取文件后缀
            self.file_ext = file_ext
            file_name = file_path.split('.')[0]  # 获取文件名称
            self.file_name = file_name

            read_funcs = {
                "csv": lambda x: pd.read_csv(x, nrows=1),
                "xls": lambda x: pd.read_excel(x, nrows=1),
                "xlsx": lambda x: pd.read_excel(x, nrows=1),
                "txt": lambda x: pd.read_table(x, nrows=1),
                "json": lambda x: pd.read_json(x, nrows=1),
                "html": lambda x: pd.read_html(x)[0].head(1)
            }

            data = read_funcs.get(file_ext)(
                file_path) if file_ext in read_funcs else None
            self.combo03_3.clear()
            self.combo03_3.addItems(["无","智能识别", *data.columns.tolist()])

    # 定义一个选择文件的方法
    def select_file(self):
        # 弹出一个文件对话框，返回文件路径
        file_path, _ = QFileDialog.getOpenFileName(self, '选择文件')
        if file_path:
            self.file_path = file_path
            self.edit02_1.setText(file_path)
            file_ext = file_path.split('.')[-1]  # 获取文件后缀
            self.file_ext = file_ext
            file_name = file_path.split('.')[0]  # 获取文件名称
            self.file_name = file_name

            read_funcs = {
                "csv": lambda x: pd.read_csv(x, nrows=1),
                "xls": lambda x: pd.read_excel(x, nrows=1),
                "xlsx": lambda x: pd.read_excel(x, nrows=1),
                "txt": lambda x: pd.read_table(x, nrows=1),
                "json": lambda x: pd.read_json(x, nrows=1),
                "html": lambda x: pd.read_html(x)[0].head(1)
            }

            data = read_funcs.get(file_ext)(
                file_path) if file_ext in read_funcs else None
            self.combo02_2.clear()
            self.combo02_2.addItems(["无", *data.columns.tolist()])
            self.combo02_3.clear()
            self.combo02_3.addItems([*data.columns.tolist()])
            
    def select_source_folder(self):
        # 弹出选择文件夹对话框，并将路径显示在输入框中
        folder_path = QFileDialog.getExistingDirectory(self, '选择原始文件夹')
        if folder_path:
            self.edit01_1.setText(folder_path)

    # 定义 open_url 槽函数，用 QDesktopServices 打开浏览器
    def open_url(self,url):
        QDesktopServices.openUrl(QUrl(url))
    def show_help(self):
        # 在这个槽函数中打开一个帮助文档或者显示一个帮助消息
        # 获取 Python 文件和 readme.txt 文件的绝对路径
        py_file = os.path.abspath(__file__)
        txt_file = os.path.join(os.path.dirname(py_file), 'readme.txt')
        # 以只读模式打开 readme.txt 文件，并读取其内容
        f = open(txt_file, 'r', encoding="UTF-8")
        text = f.read()
        f.close()
        # 创建一个消息框对象，并设置标题，图标，和文本
        msg = QMessageBox()
        msg.setWindowTitle('帮助文档')
        msg.setIcon(QMessageBox.Information)
        msg.setText(text)
        # 以模态的方式显示消息框对象
        msg.exec_()
        
    def export_excel(self):
        try:
            self.df.to_excel(self.path,na_rep="",index=False)
        except:
            self.df.to_csv(self.path,na_rep="",index=False)

        self.subWindow01_sub_window.close()
        QMessageBox.information(self, '提示', f'成功合并{self.count}个文件!')

    # 定义事件处理函数：打开子窗口
    def open_subwindow(self):
        self.subWindow01_sub_window.show()
        df = self.df.fillna("")
        nr, nc = df.shape # 获取行数和列数 
        self.table.setRowCount(nr) # 设置表格行数
        self.table.setColumnCount(nc) # 设置表格列数 
        self.table.setHorizontalHeaderLabels(df.columns)
        for r in range(nr): # 遍历每一行 
            for c in range(nc): # 遍历每一列 
                value = df.iloc[r,c] # 获取单元格值 
                item = QTableWidgetItem(str(value)) # 创建单元格对象 
                self.table.setItem(r,c,item) # 设置单元格内容
    def export_csv(self):
        try:
            to_funcs = {
                    "csv": self.df02.to_csv,
                    "xls": self.df02.to_excel,
                    "xlsx": self.df02.to_excel,
                }

            to_funcs.get(self.file_ext)(
                self.save_path) if self.file_ext in to_funcs else None
            QMessageBox.information(self, '提示', '导出文件成功!')
            # return
        except:
            QMessageBox.information(self, '提示', '导出文件失败!')

    def export_html(self):
        self.html_path = self.save_path
        try:
            self.report.to_file(self.html_path)
            QMessageBox.information(self, '提示', '导出html成功!')
            # return
        except:
            QMessageBox.information(self, '提示', '导出html失败!')

    # 定义一个槽函数，根据选项切换图形
    def switch_figure(self):
        # 获取当前选中的选项
        option = self.combo02_2_layout_subWindow02_sub_window.currentText()
        layout_subWindow02_sub_window = self.layout_subWindow02_sub_window
        fig01 = self.fig01
        fig02 = self.fig02
        # 如果选中的是条形图，那么展示fig01
        if option == "条形图":
            layout_subWindow02_sub_window.addWidget(fig01, 2, 0,1,10) # 将fig01添加到布局中的第二行第一列
            fig02.hide() # 隐藏fig02
            fig01.show() # 显示fig01
        # 如果选中的是饼状图，那么展示fig02
        elif option == "饼状图":
            layout_subWindow02_sub_window.addWidget(fig02, 2, 0,1,10) # 将fig02添加到布局中的第二行第一列
            fig01.hide() # 隐藏fig01
            fig02.show() # 显示fig02

    def show_subwindow(self):
        # 创建子窗口对象
        self.subWindow02_sub_window = QWidget()
        self.subWindow02_sub_window.setWindowTitle("同源分组")
        # self.subWindow02_sub_window.resize(900, 1200)
        # 创建子窗口中的确认按钮对象
        self.confirm_button02_1 = QPushButton("导出", self.subWindow02_sub_window)
        self.confirm_button02_1.clicked.connect(self.export_csv)

        # 创建一个网格布局，将网页视图、按钮和标签放入布局中
        self.layout_subWindow02_sub_window = QGridLayout(self.subWindow02_sub_window)
        
        # 添加图形选项
        self.label_option02_2_layout_subWindow02_sub_window = QLabel("图形类型:")
        self.combo02_2_layout_subWindow02_sub_window = QComboBox()
        self.combo02_2_layout_subWindow02_sub_window.addItem("条形图")
        self.combo02_2_layout_subWindow02_sub_window.addItem("饼状图")

        # QComboBox控件选项连接槽函数
        self.combo02_2_layout_subWindow02_sub_window.currentTextChanged.connect(lambda:self.switch_figure())
        # 设置布局中边距
        self.layout_subWindow02_sub_window.setContentsMargins(10, 10, 10, 10)
        self.layout_subWindow02_sub_window.addWidget(self.confirm_button02_1, 0, 0, 1, 1)
        self.layout_subWindow02_sub_window.addWidget(self.label_option02_2_layout_subWindow02_sub_window, 1, 0,1,1)
        self.layout_subWindow02_sub_window.addWidget(self.combo02_2_layout_subWindow02_sub_window, 1, 1,1,1)

        # 默认显示条形图
        self.switch_figure()
        # 设置窗口的布局
        self.setLayout(self.layout_subWindow02_sub_window)
        self.subWindow02_sub_window.show()

    def show_subwindow02(self):
        self.temp_file = tempfile.NamedTemporaryFile(suffix=".html", delete=False)
        self.temp_file.write(self.html.encode("utf-8"))
        self.temp_file.close()
        # 创建子窗口对象
        self.subWindow03_sub_window = QWidget()
        self.subWindow03_sub_window.setWindowTitle("数据分析")
        self.subWindow03_sub_window.resize(900, 1000)
        # 创建子窗口中的确认按钮对象
        self.confirm_button03_1 = QPushButton("导出", self.subWindow03_sub_window)
        self.confirm_button03_1.clicked.connect(self.export_html)

        # 创建一个网页视图，用于显示数据报告
        self.webview = QWebEngineView()
        # 设置默认的文本编码为UTF-8
        self.webview.settings().setDefaultTextEncoding("utf-8")

        # 创建一个网格布局，将网页视图、按钮和标签放入布局中
        self.layout_subWindow03_sub_window = QGridLayout(self.subWindow03_sub_window)
        # 设置布局中边距
        self.layout_subWindow03_sub_window.setContentsMargins(10, 10, 10, 10)
        self.layout_subWindow03_sub_window.addWidget(self.confirm_button03_1, 0, 0, 1, 1)
        self.layout_subWindow03_sub_window.addWidget(self.webview, 1, 0,1,10)
        # 设置窗口的布局
        self.setLayout(self.layout_subWindow03_sub_window)
        self.subWindow03_sub_window.show()
        # 使用QWebEngineView::load()方法加载临时文件  
        self.webview.load(QUrl.fromLocalFile(self.temp_file.name))
        self.webview.loadFinished.connect(self.delete_temp_file)
    def delete_temp_file(self):
        # 删除临时文件
        os.remove(self.temp_file.name)
    def show_view(self,result):
        if self.model == 0:
            self.df = result.df
            self.path = result.path
            self.count = result.count
            self.open_subwindow()
        elif self.model == 1:
            self.save_path = result.save_path
            self.file_ext = result.file_ext
            self.df02 = result.df02
            self.fig01 = result.fig01
            self.fig02 = result.fig02
            self.show_subwindow()
        elif self.model == 2:
            self.html = result.html
            self.save_path = result.save_path
            self.file_ext = result.file_ext
            self.report = result.report
            self.show_subwindow02()
            
    # 定义一个按钮点击函数
    def click(self, model):
        self.model = model
        data = MainData().get_data(self)
        result = MainLogic().process_data(self,data)
        self.show_view(result)

class MainData:
    def __init__(self):
        self.file_type = EXCEL_TYPE
    def get_data(self,gui):
        # 获取用户输入的参数
        if gui.model == 0:
            try:
                source_folder = Path(gui.edit01_1.text())
                if not source_folder.name:
                    QMessageBox.information(gui, '提示', '请选择原始文件夹路径!')
                    return
            except:
                QMessageBox.information(gui, '提示', '请选择原始文件夹路径!')
                return
            self.source_folder = source_folder

            try:
                is_drop_duplicates = str(gui.combo01_2.currentText())
                is_drop_duplicates = True if is_drop_duplicates=="是" else False
            except:
                is_drop_duplicates = True
            self.is_drop_duplicates = is_drop_duplicates
            # fitst_key默认None
            try:
                first_key = str(gui.edit01_3.text()).strip()
            except:
                first_key = None
            self.first_key = first_key
            # second_key默认None
            try:
                second_key = str(gui.edit01_4.text()).strip()
            except:
                second_key = None
            self.second_key = second_key
        elif gui.model == 1:

            str_fillna = False
            self.str_fillna = str_fillna
            value_fillna = False
            self.value_fillna = value_fillna
            try:
                first_class_name = str(gui.combo02_2.currentText()).strip()
                if first_class_name == "无":
                    first_class_name = False
            except:
                first_class_name = False
            self.first_class_name = first_class_name
            try:
                select_features = str(gui.combo02_3.currentText()).strip()
                if select_features == "无":
                    select_features = False
            except:
                select_features = False
            self.select_features = select_features
            try:
                n_clusters = int(gui.edit02_4.text())
            except:
                n_clusters = 2
            self.n_clusters = n_clusters
            # 如果文件路径不为空，根据文件后缀判断文件类型，并读取文件
            try:
                file_path = str(gui.edit02_1.text())
                if not file_path:
                    QMessageBox.information(gui, '提示', '请选择数据文件!')
                    return
                # file_ext = file_path.suffix # 获取文件后缀
                file_ext = file_path.split('.')[-1]  # 获取文件后缀
                self.file_ext = file_ext
                # file_name = file_path.with_suffix("")  # 获取文件名称
                file_name = file_path.split('.')[0]  # 获取文件名称
                self.file_name = file_name
                self.save_path = file_name+"_"+str(datetime.datetime.now().strftime('%Y%m%d%H%M%S'))+"."+file_ext
                # self.save_path = file_name.with_name(file_name.stem + "_group").with_suffix(file_ext)
                
                read_funcs = {
                    "csv": pd.read_csv,
                    "xls": pd.read_excel,
                    "xlsx": pd.read_excel,
                    "txt": pd.read_table,
                    "json": pd.read_json,
                    "html": lambda x: pd.read_html(x)[0]
                }
                data = read_funcs.get(file_ext)(
                    file_path) if file_ext in read_funcs else None

                # 如果数据为空，显示错误信息
                if data is not None:
                    self.data = data
                else:
                    QMessageBox.information(self, '提示', '数据为空!')
                    return
            except:
                QMessageBox.information(gui, '提示', '请选择数据文件!')
                return
            
        elif gui.model == 2:
            # str_fillna默认False
            try:
                str_fillna = str(gui.combo03_2.currentText()).strip()
            except:
                str_fillna = False
            self.str_fillna = str_fillna
            # value_fillna默认False
            try:
                if self.str_fillna == "填充固定值":
                    value_fillna = float(gui.edit03_2.text())
                    if not value_fillna:
                        value_fillna = 0
                else:
                    value_fillna = False
            except:
                value_fillna = False
            self.value_fillna = value_fillna
            try:
                timeseries_name = str(gui.combo03_3.currentText()).strip()
                if timeseries_name == "无":
                    timeseries_name = False
            except:
                timeseries_name = False
            self.timeseries_name = timeseries_name
            # 如果文件路径不为空，根据文件后缀判断文件类型，并读取文件
            try:
                file_path = str(gui.edit03_1.text())
                if not file_path:
                    QMessageBox.information(gui, '提示', '请选择数据文件!')
                    return
                # file_ext = file_path.suffix # 获取文件后缀
                file_ext = file_path.split('.')[-1]  # 获取文件后缀
                self.file_ext = file_ext
                # file_name = file_path.with_suffix("")  # 获取文件名称
                file_name = file_path.split('.')[0]  # 获取文件名称
                self.file_name = file_name
                self.save_path = file_name+"_"+str(datetime.datetime.now().strftime('%Y%m%d%H%M%S'))+".html"
                
                read_funcs = {
                    "csv": pd.read_csv,
                    "xls": pd.read_excel,
                    "xlsx": pd.read_excel,
                    "txt": pd.read_table,
                    "json": pd.read_json,
                    "html": lambda x: pd.read_html(x)[0]
                }
                data = read_funcs.get(file_ext)(
                    file_path) if file_ext in read_funcs else None

                # 如果数据为空，显示错误信息
                if data is not None:
                    self.data = data
                else:
                    QMessageBox.information(self, '提示', '数据为空!')
                    return
            except:
                QMessageBox.information(gui, '提示', '请选择数据文件!')
                return

        return self

class MainLogic:
    def __init__(self):
        pass

    # 定义一个函数，计算两个氨基酸序列的相似度，采用全局序列比对的方法
    @staticmethod
    def similarity(seq1, seq2):
        # 创建一个PairwiseAligner对象，设置全局序列比对的参数
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        aligner.match_score = 1
        aligner.mismatch_score = 0
        aligner.open_gap_score = -0.5
        aligner.extend_gap_score = -0.1
        # 使用PairwiseAligner对象进行序列比对，返回一个Alignment对象，包含所有可能的比对结果
        alignments = aligner.align(seq1, seq2)
        # 选择最高得分的比对结果
        best_score = alignments.score
        # 返回相似度，即最高得分除以较长序列的长度
        return best_score / max(len(seq1), len(seq2))
    
    # 定义一个函数，计算一行或一列的相似度矩阵，并返回一个numpy数组
    @staticmethod
    def similarity_row_or_column(i,sequences):
        # sequences = self.sequences
        # 创建一个空数组，用来存储一行或一列的相似度值
        sim_row_or_col = np.zeros(len(sequences))
        
        # 获取当前序列（即当前样本的fwr1特征）
        seq1 = sequences[i]
        
        # 遍历sequences列表中从当前序列开始的后续序列（即后续样本的fwr1特征）
        for j in range(i, len(sequences)):
            # 获取后续序列（即后续样本的fwr1特征）
            seq2 = sequences[j]
            
            # 调用similarity函数，计算当前序列和后续序列之间的相似度
            sim = MainLogic.similarity(seq1, seq2)
            
            # 将相似度存储到数组中（对称矩阵）
            sim_row_or_col[i] = sim
            sim_row_or_col[j] = sim
        
        # 返回数组
        return sim_row_or_col
    
    def similarity_sequences(self,array):
        # 使用进程池并行计算相似度矩阵，并将结果转换为numpy矩阵
        for i in range(len(array)):
            # 创建一个进程池，设置进程数为CPU核心数
            pool = Pool()
            sequences = array[i]
            similarity_matrix = np.array(pool.starmap(MainLogic.similarity_row_or_column, zip(
            range(len(sequences)), [sequences]*len(sequences))))
            # 关闭进程池，释放资源
            pool.close()
            pool.join()
            # 定义一个变量，表示要分成的组数（10组）
            n_clusters = self.n_clusters

            # 创建一个KMeans对象，设置聚类参数，如组数，随机种子等
            if len(similarity_matrix) < 10000:
                kmeans = KMeans(n_clusters=n_clusters, random_state=0,algorithm="auto",init="k-means++")
            else:
                kmeans = MiniBatchKMeans(n_clusters=n_clusters, batch_size=100,random_state=0,init="k-means++")

            # 使用KMeans对象对相似度矩阵进行聚类分析，返回一个数组，包含每个样本所属的组编号（从0到9）
            labels = kmeans.fit_predict(similarity_matrix)

            # 在df01中添加一个标签label，值为对应的组编号（将labels数组转换为一个Series对象，并赋值给df01的label列）
            self.label_name = f"label_{self.select_features}"
            self.df_dic[self.first_class[i]][self.label_name] = pd.Series(labels)
        return self.df_dic
    def process_data(self,gui,data0):
        if gui.model == 0:
            count = 0
            file_type_lower = [(".%s" % x.lower()) for x in data0.file_type]
            df_all = pd.DataFrame([])
            self.is_csv = False
            for file in data0.source_folder.iterdir():
                if file.is_file() and ((file.suffix).lower() in file_type_lower):
                    try:
                        if (file.suffix).lower() == ".csv":
                            df_file = pd.read_csv(file)
                            self.is_csv = True
                        else:
                            df_file = pd.read_excel(file)
                    except:
                        pass
                    df_file.dropna(how="all", axis="index")
                    # 获取文件的创建时间戳
                    ctime = os.path.getctime(file)
                    ctime_dt = datetime.datetime.fromtimestamp(ctime)
                    # 获取文件的修改时间戳
                    mtime = os.path.getmtime(file)
                    mtime_dt = datetime.datetime.fromtimestamp(mtime)
                    df_file["time"] = ctime_dt
                    df_file["time_m"] = mtime_dt
                    df_all = pd.concat([df_all, df_file],
                                    axis="index", ignore_index=True)
                    count += 1
            key_list = []
            first_key = data0.first_key
            if not first_key:
                first_key = df_all.columns[0]
                key_list.append(first_key)
            else:
                if first_key in df_all.columns:
                    key_list.append(first_key)
                else:
                    QMessageBox.information(gui, '提示', '去重列名错误!')
                    return
                
            second_key=data0.second_key
            if second_key:
                if (second_key in df_all.columns):
                    key_list.append(second_key)
                else:
                    QMessageBox.information(gui, '提示', '去重列名错误!')
                    return
            is_drop_duplicates = data0.is_drop_duplicates
            if is_drop_duplicates:
                df_all = df_all.sort_values(by=["time_m","time"],ascending=False).drop_duplicates(subset=key_list,keep="first").sort_index(axis="index")
            df_all.drop(columns=["time_m","time"],inplace=True)

            new_name = "merge_"+str(datetime.datetime.now().strftime('%Y%m%d%H%M%S'))
            new_path = data0.source_folder / new_name
            self.df = df_all
            if self.is_csv:
                self.path = new_path.with_suffix(".csv")
            else:
                self.path = new_path.with_suffix(".xlsx")
            self.count = count

        elif gui.model == 1:
            name = os.path.basename(gui.file_path)
            data = data0.data.copy()

            self.select_features = data0.select_features
            data.dropna(subset=[self.select_features],inplace=True)
            first_class_name = data0.first_class_name
            if first_class_name == "无":
                first_class_name = False

            self.n_clusters = data0.n_clusters
            if not first_class_name:
                first_class_name = "all_class"
                data[first_class_name] = "all"

            self.first_class=data[first_class_name].unique()
            if len(self.first_class)>1:
                self.df_dic = {name:group.reset_index() for name,group in data.groupby(first_class_name)}
            else:
                self.df_dic = {self.first_class[0]:data}
            self.array_sequences = np.array([(group[self.select_features]) for group in self.df_dic.values()],dtype=object)

            self.save_path = data0.save_path
            self.file_ext = data0.file_ext

            # 使用mmap模块创建一个内存映射文件对象，将array_sequences写入到内存中，并返回其地址和大小
            with open(data0.file_name+".dat", "wb") as f:
                dill.dump(self.array_sequences, f) # 使用dill模块的dump函数将array_sequences序列化并写入文件
            with open(data0.file_name+".dat", "r+b") as f:
                self.mm_array_sequences = mmap.mmap(f.fileno(), 0)

            # 调用similarity_sequences函数，并传入内存映射文件对象作为参数
            rev_seq=dill.load(self.mm_array_sequences)
            self.df_dic02 = self.similarity_sequences(rev_seq) # 使用dill模块的load函数将mm_array_sequences反序列化为array_sequences
            self.mm_array_sequences.close()
            f.close()
            os.remove(data0.file_name+".dat")

            self.df02 = pd.concat(self.df_dic02.values(),axis=0)
            try:
                df_g=self.df02.groupby([first_class_name,self.label_name]).size().reset_index(name="count")
                df_g_p=df_g.pivot(index=first_class_name,columns=self.label_name,values="count")
            except:
                QMessageBox.information(gui, '提示', '绘图错误!')
                return

            fig01, ax = plt.subplots()
            df_g_p.plot.bar(color=[matplotlib.colors.to_rgb((random.uniform(0,1),random.
            uniform(0,1),random.uniform(0,1))) for _ in range(self.n_clusters)],title="counts", ax=ax)
            # 添加文本标签
            for p in ax.patches: # 遍历每个柱形对象
                x = p.get_x() + p.get_width() / 2 # 获取x坐标
                y = p.get_y() + p.get_height() # 获取y坐标
                h = p.get_height() # 获取高度
                ax.text(x, y, h, ha='center', va='bottom', fontsize=10) # 在上方添加文本标签

            # 获取first_class_name的唯一值列表
            fc_values = np.array(df_g_p.index.unique())

            # 创建一个n行1列的子图，其中n是first_class_name的唯一值的个数+1
            fig02, axes = plt.subplots(1,len(fc_values))

            # 遍历每个子图和每个first_class_name的值
            if len(fc_values) > 1:
                for ax, fc in zip(axes, fc_values):
                    # 根据该值筛选出对应的数据
                    data = df_g_p.loc[fc]
                    labels = data.index # 获取每个类别的标签
                    colors = [matplotlib.colors.to_rgb((random.uniform(0,1),random.
                    uniform(0,1),random.uniform(0,1))) for _ in range(self.n_clusters)] # 随机生成颜色列表
                    # 在子图中绘制饼状图
                    ax.pie(data, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False) # 绘制饼状图
                    # 设置子图的标题
                    ax.set_title(f"{first_class_name} = {fc}")
            else:
                # 根据该值筛选出对应的数据
                data = df_g_p.loc["all"]
                labels = data.index # 获取每个类别的标签
                colors = [matplotlib.colors.to_rgb((random.uniform(0,1),random.
                uniform(0,1),random.uniform(0,1))) for _ in range(self.n_clusters)] # 随机生成颜色列表
                # 在子图中绘制饼状图
                axes.pie(data, labels=labels, colors=colors, autopct='%1.1f%%', shadow=False) # 绘制饼状图
                # 设置子图的标题
                # axes.set_title(f"{first_class_name} = {fc_values[0]}")

            # 调整子图间距
            plt.tight_layout() # 自动调整子图间距

            self.fig01 = FigureCanvasQTAgg(fig01)
            self.fig02 = FigureCanvasQTAgg(fig02)

        elif gui.model == 2:
            name = os.path.splitext(os.path.basename(data0.save_path))[0]
            data = data0.data.copy()
            self.save_path = data0.save_path
            self.file_ext = data0.file_ext

            str_fillna = data0.str_fillna
            value_fillna = data0.value_fillna
            if str_fillna == "去除缺失值":
                data.dropna(axis=0, how="any", inplace=True)
            elif str_fillna == '填充固定值':
                data.fillna(value=value_fillna, inplace=True)
            elif str_fillna == '填充均值':
                for col in data.columns:
                    if data[col].dtype in ['int', 'float']:
                        data[col].fillna(data[col].mean(), inplace=True)
            elif str_fillna == '填充中位数':
                for col in data.columns:
                    if data[col].dtype in ['int', 'float']:
                        data[col].fillna(
                            data[col].median(), inplace=True)
            elif str_fillna == '填充前值':
                data.fillna(method="ffill", inplace=True)
            elif str_fillna == '填充后值':
                data.fillna(method="bfill", inplace=True)
            tsmode = False
            sortby = None  # 给sortby一个初始值
            timeseries_name = data0.timeseries_name
            try:
                if timeseries_name:
                    tsmode = True
                    sortby = timeseries_name
                    if sortby == "智能识别":
                        sortby = None
                        for col in data.columns:
                            try:
                                data[col] = pd.to_datetime(data[col])
                                if data[col].dtype == "datetime64[ns]":
                                    sortby = col
                                    break
                            except:
                                pass
                    else:
                        data[sortby] = pd.to_datetime(data[sortby])
            except:
                tsmode = False
            minimal = True if data.shape[0] > 100000 or data.shape[1] > 100 else False
            if sortby is not None:
                self.report = ydata_profiling.ProfileReport(data, title=str(
                    name), tsmode=tsmode, minimal=minimal, sortby=sortby)
            else:
                self.report = ydata_profiling.ProfileReport(
                    data, title=str(name), tsmode=tsmode, minimal=minimal)
            self.html = self.report.to_html()  # 将数据报告转换为HTML格式的字符串
        return self

if __name__ == '__main__':  # 主程序入口
    # 调用multiprocessing.freeze_support()函数
    multiprocessing.freeze_support()
    # 创建应用对象
    app = QApplication(sys.argv)
    # 创建窗口对象并显示 
    window = MainView()
    window.show()
    # 进入事件循环
    sys.exit(app.exec_())