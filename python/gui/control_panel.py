from typing import Any

from PyQt6.QtCore import QLocale, QRegularExpression, Qt, pyqtSignal
from PyQt6.QtGui import QRegularExpressionValidator
from PyQt6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QSlider,
    QSpinBox,
    QVBoxLayout,
)

from gui.setup_schema import (
    INT_PARAMS,
    PARAM_LABELS,
    PARAM_RANGES,
    PARAMS_BY_SETUP,
    SETUP_DEFAULTS,
    TEXT_PARAMS,
    decimals_for_param,
    step_for_param,
)
from gui.ui_styles import (
    StyledButton,
    checkbox_style,
    combo_style,
    panel_style,
    section_label_style,
    slider_style,
    small_label_style,
    spin_style,
)


class ControlPanel(QFrame):
    configuration_applied = pyqtSignal(str, dict)
    view_options_changed = pyqtSignal()

    def __init__(self):
        super().__init__()

        self.param_widgets: dict[str, QSpinBox | QDoubleSpinBox | QLineEdit] = {}

        self.setFixedWidth(310)
        self.setStyleSheet(panel_style())

        self._setup_ui()

    # ------------------------------------------------------------------
    # Public accessors used by MainWindow
    # ------------------------------------------------------------------

    def selected_field(self) -> str:
        return self.field_selector.currentText()

    def steps_per_tick(self) -> int:
        return self.steps_slider.value()

    def arrows_enabled(self) -> bool:
        return self.chk_arrows.isChecked()

    def uniform_arrows_enabled(self) -> bool:
        return self.chk_uniform_arrows.isChecked()

    def arrow_density(self) -> int:
        return self.arrow_density_slider.value()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _setup_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(10)

        layout.addWidget(self._section_label("SETUP"))

        self.setup_selector = QComboBox()
        self.setup_selector.addItems(list(SETUP_DEFAULTS.keys()))
        self.setup_selector.setStyleSheet(combo_style())
        self.setup_selector.currentTextChanged.connect(self._on_setup_changed)

        layout.addWidget(self.setup_selector)

        layout.addWidget(self._section_label("PARAMETERS"))

        self.param_frame = QFrame()
        self.param_frame.setStyleSheet("QFrame { border: none; }")

        self.param_layout = QFormLayout(self.param_frame)
        self.param_layout.setContentsMargins(0, 0, 0, 0)
        self.param_layout.setSpacing(8)

        initial_setup = self.setup_selector.currentText()
        self._rebuild_param_widgets(SETUP_DEFAULTS[initial_setup])

        layout.addWidget(self.param_frame)

        self.btn_apply = StyledButton("Apply setup", accent="#8be9fd")
        self.btn_apply.clicked.connect(self._emit_configuration)
        layout.addWidget(self.btn_apply)

        layout.addWidget(self._separator())
        layout.addWidget(self._section_label("VIEW"))

        layout.addLayout(self._build_field_row())

        self.chk_arrows = QCheckBox("Velocity vectors")
        self.chk_arrows.setChecked(True)
        self.chk_arrows.setStyleSheet(checkbox_style())
        self.chk_arrows.stateChanged.connect(self.view_options_changed.emit)
        layout.addWidget(self.chk_arrows)

        self.chk_uniform_arrows = QCheckBox("Uniform arrow length")
        self.chk_uniform_arrows.setChecked(False)
        self.chk_uniform_arrows.setStyleSheet(checkbox_style())
        self.chk_uniform_arrows.stateChanged.connect(self.view_options_changed.emit)
        layout.addWidget(self.chk_uniform_arrows)

        layout.addLayout(self._build_arrow_density_row())

        layout.addWidget(self._separator())
        layout.addWidget(self._section_label("SIMULATION"))

        layout.addLayout(self._build_steps_row())

        layout.addStretch()

    def _build_field_row(self) -> QHBoxLayout:
        row = QHBoxLayout()

        row.addWidget(self._small_label("Field:"))

        self.field_selector = QComboBox()
        self.field_selector.addItems(["Pressure", "Speed"])
        self.field_selector.setStyleSheet(combo_style())
        self.field_selector.currentTextChanged.connect(self.view_options_changed.emit)

        row.addWidget(self.field_selector)

        return row

    def _build_arrow_density_row(self) -> QHBoxLayout:
        row = QHBoxLayout()

        row.addWidget(self._small_label("Arrow density:"))

        self.arrow_density_slider = QSlider(Qt.Orientation.Horizontal)
        self.arrow_density_slider.setRange(8, 80)
        self.arrow_density_slider.setValue(30)
        self.arrow_density_slider.setStyleSheet(slider_style())
        self.arrow_density_slider.valueChanged.connect(
            self._on_arrow_density_changed
        )

        self.arrow_density_label = QLabel("30")
        self.arrow_density_label.setStyleSheet(small_label_style())

        row.addWidget(self.arrow_density_slider)
        row.addWidget(self.arrow_density_label)

        return row

    def _build_steps_row(self) -> QHBoxLayout:
        row = QHBoxLayout()

        row.addWidget(self._small_label("Steps/tick:"))

        self.steps_slider = QSlider(Qt.Orientation.Horizontal)
        self.steps_slider.setRange(1, 10)
        self.steps_slider.setValue(1)
        self.steps_slider.setStyleSheet(slider_style())

        self.steps_label = QLabel("1")
        self.steps_label.setStyleSheet(small_label_style())

        self.steps_slider.valueChanged.connect(self._on_steps_changed)

        row.addWidget(self.steps_slider)
        row.addWidget(self.steps_label)

        return row

    # ------------------------------------------------------------------
    # Parameter widgets
    # ------------------------------------------------------------------

    def _rebuild_param_widgets(self, params: dict[str, Any]) -> None:
        while self.param_layout.rowCount() > 0:
            self.param_layout.removeRow(0)

        self.param_widgets.clear()

        setup_name = self.setup_selector.currentText()
        param_names = PARAMS_BY_SETUP[setup_name]

        for name in param_names:
            label = self._small_label(PARAM_LABELS.get(name, name))
            value = params[name]

            if name in TEXT_PARAMS:
                widget = self._make_text_param(name, str(value))

            elif name in INT_PARAMS:
                widget = self._make_int_param(name, int(value))

            else:
                widget = self._make_float_param(name, float(value))

            self.param_widgets[name] = widget
            self.param_layout.addRow(label, widget)

    def _make_text_param(self, name: str, value: str) -> QLineEdit:
        widget = QLineEdit()
        widget.setText(value)
        widget.setStyleSheet(spin_style())

        if name == "naca_code":
            widget.setMaxLength(4)
            widget.setPlaceholderText("2412")

            validator = QRegularExpressionValidator(
                QRegularExpression(r"\d{0,4}")
            )
            widget.setValidator(validator)

        return widget

    def _make_int_param(self, name: str, value: int) -> QSpinBox:
        widget = QSpinBox()
        widget.setLocale(
            QLocale(
                QLocale.Language.English,
                QLocale.Country.UnitedStates,
            )
        )

        widget.setRange(1, 10_000)
        widget.setSingleStep(1)
        widget.setValue(value)
        widget.setStyleSheet(spin_style())

        return widget

    def _make_float_param(self, name: str, value: float) -> QDoubleSpinBox:
        widget = QDoubleSpinBox()
        widget.setLocale(
            QLocale(
                QLocale.Language.English,
                QLocale.Country.UnitedStates,
            )
        )

        lo, hi = PARAM_RANGES.get(name, (1e-9, 1e9))
        widget.setRange(lo, hi)
        widget.setDecimals(decimals_for_param(name))
        widget.setSingleStep(step_for_param(name))
        widget.setValue(value)
        widget.setStyleSheet(spin_style())

        return widget

    def _read_params(self) -> dict[str, Any]:
        params: dict[str, Any] = {}

        for name, widget in self.param_widgets.items():
            if isinstance(widget, QLineEdit):
                params[name] = widget.text().strip()

            elif name in INT_PARAMS:
                params[name] = int(widget.value())

            else:
                params[name] = float(widget.value())

        return params

    # ------------------------------------------------------------------
    # Events
    # ------------------------------------------------------------------

    def _on_setup_changed(self, setup_name: str) -> None:
        self._rebuild_param_widgets(SETUP_DEFAULTS[setup_name])

    def _emit_configuration(self) -> None:
        setup_name = self.setup_selector.currentText()
        params = self._read_params()

        self.configuration_applied.emit(setup_name, params)

    def _on_arrow_density_changed(self, value: int) -> None:
        self.arrow_density_label.setText(str(value))
        self.view_options_changed.emit()

    def _on_steps_changed(self, value: int) -> None:
        self.steps_label.setText(str(value))

    # ------------------------------------------------------------------
    # Small widget helpers
    # ------------------------------------------------------------------

    def _section_label(self, text: str) -> QLabel:
        label = QLabel(text)
        label.setStyleSheet(section_label_style())
        return label

    def _small_label(self, text: str) -> QLabel:
        label = QLabel(text)
        label.setStyleSheet(small_label_style())
        return label

    def _separator(self) -> QFrame:
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet("color: #2a2f4a;")
        return sep