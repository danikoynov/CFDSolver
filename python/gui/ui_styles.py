from PyQt6.QtWidgets import QPushButton


def main_window_style() -> str:
    return "QMainWindow { background: #0e0f18; }"


def panel_style() -> str:
    return (
        "QFrame {"
        "  background: #13141f;"
        "  border: 1px solid #2a2f4a;"
        "  border-radius: 6px;"
        "}"
    )


def title_label_style() -> str:
    return (
        "color: #6272a4;"
        "font-family: 'JetBrains Mono', 'Consolas', monospace;"
        "font-size: 13px;"
        "font-weight: bold;"
        "letter-spacing: 0.18em;"
    )


def section_label_style() -> str:
    return (
        "color: #6272a4;"
        "font-family: 'JetBrains Mono', 'Consolas', monospace;"
        "font-size: 11px;"
        "font-weight: bold;"
        "letter-spacing: 0.16em;"
        "border: none;"
    )


def small_label_style() -> str:
    return (
        "color: #8899bb;"
        "font-family: 'Consolas', monospace;"
        "font-size: 11px;"
        "border: none;"
    )


def muted_label_style() -> str:
    return (
        "color: #555970;"
        "font-family: 'JetBrains Mono', 'Consolas', monospace;"
        "font-size: 11px;"
    )


def iteration_label_style() -> str:
    return (
        "color: #50fa7b;"
        "font-family: 'JetBrains Mono', 'Consolas', monospace;"
        "font-size: 15px;"
        "font-weight: bold;"
    )


def combo_style() -> str:
    return (
        "QComboBox {"
        "  background: #1e2030;"
        "  color: #c8d0e8;"
        "  border: 1px solid #3a3f58;"
        "  border-radius: 4px;"
        "  padding: 3px 8px;"
        "  font-family: 'Consolas', monospace;"
        "  font-size: 11px;"
        "}"
        "QComboBox::drop-down {"
        "  border: none;"
        "}"
        "QComboBox QAbstractItemView {"
        "  background: #1e2030;"
        "  color: #c8d0e8;"
        "  selection-background-color: #2a2f4a;"
        "}"
    )


def spin_style() -> str:
    return (
        "QSpinBox, QDoubleSpinBox {"
        "  background: #1e2030;"
        "  color: #c8d0e8;"
        "  border: 1px solid #3a3f58;"
        "  border-radius: 4px;"
        "  padding: 2px 4px;"
        "  font-family: 'Consolas', monospace;"
        "  font-size: 11px;"
        "}"
        "QSpinBox::up-button, QSpinBox::down-button,"
        "QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {"
        "  width: 14px;"
        "  background: #2a2f4a;"
        "  border: none;"
        "}"
    )


def checkbox_style() -> str:
    return (
        "QCheckBox {"
        "  color: #8899bb;"
        "  font-family: 'Consolas', monospace;"
        "  font-size: 11px;"
        "  spacing: 5px;"
        "}"
        "QCheckBox::indicator {"
        "  width: 13px;"
        "  height: 13px;"
        "  border: 1px solid #3a3f58;"
        "  border-radius: 2px;"
        "  background: #1e2030;"
        "}"
        "QCheckBox::indicator:checked {"
        "  background: #6272a4;"
        "}"
    )


def slider_style() -> str:
    return (
        "QSlider::groove:horizontal {"
        "  height: 4px;"
        "  background: #2a2f4a;"
        "  border-radius: 2px;"
        "}"
        "QSlider::handle:horizontal {"
        "  width: 12px;"
        "  height: 12px;"
        "  margin: -4px 0;"
        "  background: #6272a4;"
        "  border-radius: 6px;"
        "}"
        "QSlider::sub-page:horizontal {"
        "  background: #6272a4;"
        "  border-radius: 2px;"
        "}"
    )


def button_style(accent: str | None = None) -> str:
    base = (
        "QPushButton {"
        "  background: #1e2030;"
        "  color: #c8d0e8;"
        "  border: 1px solid #3a3f58;"
        "  border-radius: 4px;"
        "  padding: 4px 14px;"
        "  font-family: 'JetBrains Mono', 'Consolas', monospace;"
        "  font-size: 12px;"
        "}"
        "QPushButton:hover {"
        "  background: #2a2f4a;"
        "  border-color: #6272a4;"
        "}"
        "QPushButton:pressed {"
        "  background: #14162a;"
        "}"
        "QPushButton:disabled {"
        "  color: #555970;"
        "  border-color: #2a2f3a;"
        "}"
    )

    if accent is None:
        return base

    return base + (
        "QPushButton {"
        f"  border-color: {accent};"
        f"  color: {accent};"
        "}"
    )


class StyledButton(QPushButton):
    def __init__(self, text: str, accent: str | None = None):
        super().__init__(text)
        self.setStyleSheet(button_style(accent))