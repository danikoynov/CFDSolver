import sys
from pathlib import Path

here = Path(__file__).resolve()
python_dir = here.parent.parent
sys.path.insert(0, str(python_dir))

from PyQt6.QtWidgets import QApplication
from gui.main_window import MainWindow


def main():
    app = QApplication(sys.argv)

    window = MainWindow()
    window.resize(1280, 720)
    window.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()