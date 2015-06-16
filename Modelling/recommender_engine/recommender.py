# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'recommender.ui'
#
# Created: Wed May 20 14:28:05 2015
#      by: pyside-uic 0.2.15 running on PySide 1.2.2
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_recommender(object):
    def setupUi(self, recommender):
        recommender.setObjectName("recommender")
        recommender.resize(393, 622)
        self.centralWidget = QtGui.QWidget(recommender)
        self.centralWidget.setObjectName("centralWidget")
        self.label = QtGui.QLabel(self.centralWidget)
        self.label.setGeometry(QtCore.QRect(20, 10, 111, 16))
        self.label.setObjectName("label")
        self.Library = QtGui.QListWidget(self.centralWidget)
        self.Library.setGeometry(QtCore.QRect(20, 31, 361, 151))
        self.Library.setObjectName("Library")
        self.Recommendations = QtGui.QListWidget(self.centralWidget)
        self.Recommendations.setGeometry(QtCore.QRect(20, 210, 291, 361))
        self.Recommendations.setObjectName("Recommendations")
        self.Score = QtGui.QListWidget(self.centralWidget)
        self.Score.setGeometry(QtCore.QRect(330, 210, 51, 361))
        self.Score.setObjectName("Score")
        self.label_2 = QtGui.QLabel(self.centralWidget)
        self.label_2.setGeometry(QtCore.QRect(20, 190, 101, 16))
        self.label_2.setObjectName("label_2")
        self.label_3 = QtGui.QLabel(self.centralWidget)
        self.label_3.setGeometry(QtCore.QRect(330, 190, 47, 13))
        self.label_3.setObjectName("label_3")
        recommender.setCentralWidget(self.centralWidget)
        self.menuBar = QtGui.QMenuBar(recommender)
        self.menuBar.setGeometry(QtCore.QRect(0, 0, 393, 21))
        self.menuBar.setObjectName("menuBar")
        recommender.setMenuBar(self.menuBar)
        self.mainToolBar = QtGui.QToolBar(recommender)
        self.mainToolBar.setObjectName("mainToolBar")
        recommender.addToolBar(QtCore.Qt.TopToolBarArea, self.mainToolBar)
        self.statusBar = QtGui.QStatusBar(recommender)
        self.statusBar.setObjectName("statusBar")
        recommender.setStatusBar(self.statusBar)

        self.retranslateUi(recommender)
        QtCore.QMetaObject.connectSlotsByName(recommender)

    def retranslateUi(self, recommender):
        recommender.setWindowTitle(QtGui.QApplication.translate("recommender", "Bacteria Recommender", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("recommender", "Library", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("recommender", "Recommendations", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("recommender", "Score", None, QtGui.QApplication.UnicodeUTF8))

