import numpy as np
import pyqtgraph as pg

from PyQt6.QtCore import QTimer
from PyQt6.QtGui import QColor, QPalette
from PyQt6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from gui.arrow_renderer import ArrowRenderer
from gui.control_panel import ControlPanel
from gui.field_renderer import FieldRenderer
from gui.simulation_controller import SimulationController
from gui.ui_styles import (
    StyledButton,
    iteration_label_style,
    main_window_style,
    muted_label_style,
    title_label_style,
)


class MainWindow(QMainWindow):
    _ARROW_EVERY = 3

    def __init__(self):
        super().__init__()

        self.setWindowTitle("FlowGrid · CFD Dashboard")
        self.setMinimumSize(950, 560)

        self.sim = SimulationController()

        self.control_panel: ControlPanel
        self.field_renderer: FieldRenderer
        self.arrow_renderer: ArrowRenderer

        self._view_counter = 0

        self.timer = QTimer()
        self.timer.timeout.connect(self._on_tick)

        self._setup_palette()
        self._setup_ui()
        self._connect_signals()

        self.update_view(force_arrows=True)

    # ------------------------------------------------------------------
    # UI setup
    # ------------------------------------------------------------------

    def _setup_palette(self) -> None:
        palette = QPalette()
        palette.setColor(QPalette.ColorRole.Window, QColor("#0e0f18"))
        palette.setColor(QPalette.ColorRole.WindowText, QColor("#c8d0e8"))
        palette.setColor(QPalette.ColorRole.Base, QColor("#13141f"))
        palette.setColor(QPalette.ColorRole.Text, QColor("#c8d0e8"))

        self.setPalette(palette)
        self.setStyleSheet(main_window_style())

    def _setup_ui(self) -> None:
        central = QWidget()
        root = QVBoxLayout(central)
        root.setContentsMargins(10, 8, 10, 8)
        root.setSpacing(8)

        root.addLayout(self._build_top_bar())

        body = QHBoxLayout()
        body.setSpacing(10)

        self.control_panel = ControlPanel()
        plot_panel = self._build_plot_panel()

        body.addWidget(self.control_panel)
        body.addLayout(plot_panel, stretch=1)

        root.addLayout(body, stretch=1)
        root.addLayout(self._build_footer())

        self.setCentralWidget(central)

    def _build_top_bar(self) -> QHBoxLayout:
        bar = QHBoxLayout()
        bar.setSpacing(8)

        title = QLabel("FLOWGRID")
        title.setStyleSheet(title_label_style())

        self.btn_start = StyledButton("▶ Run", accent="#50fa7b")
        self.btn_pause = StyledButton("⏸ Pause", accent="#ffb86c")
        self.btn_reset = StyledButton("↺ Reset", accent="#ff5555")

        self.btn_start.clicked.connect(self.start_simulation)
        self.btn_pause.clicked.connect(self.pause_simulation)
        self.btn_reset.clicked.connect(self.reset_simulation)

        self.btn_pause.setEnabled(False)

        bar.addWidget(title)
        bar.addSpacing(8)
        bar.addWidget(self.btn_start)
        bar.addWidget(self.btn_pause)
        bar.addWidget(self.btn_reset)
        bar.addStretch()

        return bar

    def _build_plot_panel(self) -> QHBoxLayout:
        pg.setConfigOptions(antialias=False, useOpenGL=True)

        self.plot = pg.PlotWidget()
        self.plot.setAspectLocked(True)
        self.plot.hideAxis("left")
        self.plot.hideAxis("bottom")
        self.plot.setBackground("#0a0b14")
        self.plot.setMouseEnabled(x=True, y=True)
        self.plot.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )

        self.field_renderer = FieldRenderer(self.plot)
        self.arrow_renderer = ArrowRenderer(self.plot)

        plot_row = QHBoxLayout()
        plot_row.setSpacing(6)
        plot_row.addWidget(self.plot, stretch=1)
        plot_row.addLayout(self.field_renderer.build_colorbar())

        self._reset_plot_range()

        return plot_row

    def _build_footer(self) -> QHBoxLayout:
        footer = QHBoxLayout()
        footer.setContentsMargins(4, 2, 4, 2)

        iter_prefix = QLabel("ITERATION")
        iter_prefix.setStyleSheet(
            "color:#3a3f58;"
            "font-family:'JetBrains Mono','Consolas',monospace;"
            "font-size:11px;"
            "letter-spacing:0.12em;"
        )

        self.iter_label = QLabel("0")
        self.iter_label.setStyleSheet(iteration_label_style())

        self.range_label = QLabel("")
        self.range_label.setStyleSheet(muted_label_style())

        footer.addWidget(iter_prefix)
        footer.addSpacing(6)
        footer.addWidget(self.iter_label)
        footer.addStretch()
        footer.addWidget(self.range_label)

        return footer

    def _connect_signals(self) -> None:
        self.control_panel.configuration_applied.connect(self.apply_configuration)
        self.control_panel.view_options_changed.connect(self._force_view_update)

    # ------------------------------------------------------------------
    # Simulation controls
    # ------------------------------------------------------------------

    def start_simulation(self) -> None:
        self.btn_start.setEnabled(False)
        self.btn_pause.setEnabled(True)
        self.timer.start(33)

    def pause_simulation(self) -> None:
        self.timer.stop()
        self.btn_start.setEnabled(True)
        self.btn_pause.setEnabled(False)

    def reset_simulation(self) -> None:
        self.timer.stop()
        self.btn_start.setEnabled(True)
        self.btn_pause.setEnabled(False)

        try:
            self.sim.reset()
            self._after_simulation_recreated()
        except Exception as exc:
            self._show_error("Reset failed", exc)

    def apply_configuration(self, setup_name: str, params: dict) -> None:
        was_running = self.timer.isActive()
        self.timer.stop()

        try:
            self.sim.configure(setup_name, params)
            self._after_simulation_recreated()

        except Exception as exc:
            self._show_error("Configuration failed", exc)

            self.btn_start.setEnabled(True)
            self.btn_pause.setEnabled(False)
            return

        if was_running:
            self.btn_start.setEnabled(False)
            self.btn_pause.setEnabled(True)
            self.timer.start(33)

    def _on_tick(self) -> None:
        self.sim.step(steps=self.control_panel.steps_per_tick())
        self.update_view()

    def _after_simulation_recreated(self) -> None:
        self._view_counter = 0
        self.arrow_renderer.clear()
        self._reset_plot_range()
        self.update_view(force_arrows=True)

    def _show_error(self, title: str, exc: Exception) -> None:
        QMessageBox.critical(
            self,
            title,
            str(exc),
        )

    # ------------------------------------------------------------------
    # Rendering
    # ------------------------------------------------------------------

    def update_view(self, force_arrows: bool = False) -> None:
        solid = self.sim.solid_mask()
        field_name = self.control_panel.selected_field()
        u = v = None

        if field_name == "Pressure":
            field = self.sim.pressure_field()
        elif field_name == "Speed":
            u, v = self.sim.velocity_fields()
            field = np.hypot(u, v)
        else:
            raise ValueError(f"Unknown field: {field_name}")

        vmin, vmax = self.field_renderer.render(field, solid)

        self._view_counter += 1

        should_draw_arrows = self.control_panel.arrows_enabled() and (
            force_arrows or self._view_counter % self._ARROW_EVERY == 0
        )

        if should_draw_arrows:
            if u is None:
                u, v = self.sim.velocity_fields()
            self.arrow_renderer.draw(
                u, v, solid,
                density=self.control_panel.arrow_density(),
                uniform_length=self.control_panel.uniform_arrows_enabled(),
            )
        elif not self.control_panel.arrows_enabled():
            self.arrow_renderer.clear()

        self.iter_label.setText(f"{self.sim.step_count:,}")
        self.range_label.setText(f"{field_name}  {vmin:.4f} … {vmax:.4f}")

    def _force_view_update(self) -> None:
        self.update_view(force_arrows=True)

    def _reset_plot_range(self) -> None:
        self.plot.setXRange(0, self.sim.grid_w, padding=0)
        self.plot.setYRange(0, self.sim.grid_h, padding=0)