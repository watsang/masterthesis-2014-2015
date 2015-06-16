from PySide.QtCore import *
from PySide.QtGui import *
import sys
from PySide import QtCore, QtGui

import recommender
import pandas as pd

# Define the libraries where you want to look up your recommendations and scores
sorted_names = pd.read_csv("sorted_names.csv", header=0)
scores = pd.read_csv("scores.csv", header=0)
scores = scores.sort(axis=1)
sorted_names = sorted_names.sort(axis=1)
names = sorted_names.columns.values

class MainDialog(QMainWindow, recommender.Ui_recommender):
    
    def __init__(self, parent = None):
        super(MainDialog, self).__init__(parent)
        self.setupUi(self) 
        self.Library.addItems(names)
        
        QtCore.QObject.connect(self.Library, QtCore.SIGNAL("itemClicked(QListWidgetItem*)"), self.updateUi)
        QtCore.QObject.connect(self.Library, QtCore.SIGNAL("itemSelectionChanged()"), self.updateUi)        
        
    def updateUi(self):      
        self.Recommendations.clear()
        self.Score.clear()
        selected = self.Library.currentItem().text()
        self.Recommendations.addItems(sorted_names[selected])
        temp_score = [str(round(i,3)) for i in scores[selected]]
        self.Score.addItems(temp_score)
        
app = QApplication.instance()
form = MainDialog()
form.show()
app.exec_()