import numpy as np
import pyqtgraph as pg

from PyQt6.QtWidgets import QLabel, QSizePolicy, QVBoxLayout


_CMAP = pg.ColorMap(
    pos=np.array([0.0, 1.0]),
    color=np.array(
        [[0, 0, 255, 255], [255, 0, 0, 255]],
        dtype=np.uint8,
    ),
)

LUT_BLUE_RED = _CMAP.getLookupTable(nPts=256, alpha=False)
SOLID_COLOR = np.array([120, 120, 120], dtype=np.uint8)


def _normalise_field(
    field: np.ndarray,
    solid: np.ndarray,
) -> tuple[np.ndarray, float, float]:
    fluid_vals = field[~solid]

    if fluid_vals.size == 0:
        vmin = 0.0
        vmax = 1.0
    else:
        vmin = float(fluid_vals.min())
        vmax = float(fluid_vals.max())

    span = vmax - vmin

    if span < 1e-12:
        normed = np.full_like(field, 128.0)
    else:
        normed = (field - vmin) / span * 255.0

    normed = np.clip(normed, 0, 255).astype(np.uint8)

    return normed, vmin, vmax


def _apply_colormap(
    field: np.ndarray,
    solid: np.ndarray,
) -> tuple[np.ndarray, float, float]:
    normed, vmin, vmax = _normalise_field(field, solid)

    rgb = LUT_BLUE_RED[normed]
    rgb = rgb.copy()
    rgb[solid] = SOLID_COLOR

    return rgb, vmin, vmax


class FieldRenderer:
    def __init__(self, plot: pg.PlotWidget):
        self.plot = plot

        self.image_item = pg.ImageItem()
        self.image_item.setOpts(axisOrder="row-major", useRGBA=True)

        self.plot.addItem(self.image_item)

        self.cb_min_label: QLabel | None = None
        self.cb_max_label: QLabel | None = None

    def render(
        self,
        field: np.ndarray,
        solid: np.ndarray,
    ) -> tuple[float, float]:
        rgb, vmin, vmax = _apply_colormap(field, solid)

        self.image_item.setImage(
            rgb,
            autoLevels=False,
            levels=[0, 255],
        )

        self.update_colorbar(vmin, vmax)

        return vmin, vmax

    def build_colorbar(self) -> QVBoxLayout:
        layout = QVBoxLayout()
        layout.setSpacing(2)
        layout.setContentsMargins(0, 0, 0, 0)

        mono = "font-family:'JetBrains Mono','Consolas',monospace; font-size:10px;"

        self.cb_max_label = QLabel("1.000")
        self.cb_max_label.setStyleSheet(
            f"color:#ff4444; {mono} qproperty-alignment: AlignHCenter;"
        )

        strip = QLabel()
        strip.setFixedWidth(18)
        strip.setSizePolicy(
            QSizePolicy.Policy.Fixed,
            QSizePolicy.Policy.Expanding,
        )
        strip.setStyleSheet(
            "background: qlineargradient("
            "  x1:0, y1:0, x2:0, y2:1,"
            "  stop:0 rgb(255,0,0),"
            "  stop:1 rgb(0,0,255)"
            ");"
            "border: 1px solid #2a2f4a;"
            "border-radius: 2px;"
        )

        self.cb_min_label = QLabel("0.000")
        self.cb_min_label.setStyleSheet(
            f"color:#4444ff; {mono} qproperty-alignment: AlignHCenter;"
        )

        layout.addWidget(self.cb_max_label)
        layout.addWidget(strip, stretch=1)
        layout.addWidget(self.cb_min_label)

        return layout

    def update_colorbar(self, vmin: float, vmax: float) -> None:
        if self.cb_max_label is not None:
            self.cb_max_label.setText(f"{vmax:.3f}")

        if self.cb_min_label is not None:
            self.cb_min_label.setText(f"{vmin:.3f}")
            