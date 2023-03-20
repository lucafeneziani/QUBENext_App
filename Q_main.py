from GUI_struct import QApp
from PyQt5.QtWidgets import QApplication
import sys

app = QApplication(sys.argv)
win = QApp()

win.show()
sys.exit(app.exec_())