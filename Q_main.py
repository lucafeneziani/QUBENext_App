from GUI_struct import QApp
from PyQt5.QtWidgets import QApplication
#from PyQt5.QtGui import QPixmap
import sys


app = QApplication(sys.argv)
win = QApp()

win.show()
sys.exit(app.exec_())