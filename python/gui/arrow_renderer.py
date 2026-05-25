import numpy as np
import pyqtgraph as pg


class ArrowRenderer:
    def __init__(self, plot: pg.PlotWidget):
        self.plot = plot

        self._shaft = pg.PlotDataItem(
            [], [], pen=pg.mkPen((180, 190, 220, 200), width=1), connect="finite"
        )
        self._head = pg.PlotDataItem(
            [], [], pen=pg.mkPen((180, 190, 220, 140), width=1), connect="finite"
        )
        self.plot.addItem(self._shaft)
        self.plot.addItem(self._head)

    def clear(self) -> None:
        self._shaft.setData([], [])
        self._head.setData([], [])

    def draw(
        self,
        u: np.ndarray,
        v: np.ndarray,
        solid: np.ndarray,
        density: int,
        uniform_length: bool,
    ) -> None:
        h, w = u.shape
        stride = max(1, max(w, h) // max(1, density))

        max_speed = float(np.hypot(u, v).max())
        if max_speed < 1e-12:
            self.clear()
            return

        js = np.arange(0, h, stride)
        is_ = np.arange(0, w, stride)
        jj, ii = np.meshgrid(js, is_, indexing="ij")

        mask = ~solid[jj, ii]
        jj, ii = jj[mask], ii[mask]

        ui = u[jj, ii]
        vi = v[jj, ii]
        speed = np.hypot(ui, vi)

        valid = speed > 1e-12
        jj, ii, ui, vi, speed = jj[valid], ii[valid], ui[valid], vi[valid], speed[valid]

        if len(ii) == 0:
            self.clear()
            return

        arrow_len = 0.42 * stride
        scale = speed if uniform_length else max_speed
        dx = ui / scale * arrow_len
        dy = vi / scale * arrow_len

        x0, y0 = ii + 0.5, jj + 0.5
        x1, y1 = x0 + dx, y0 + dy

        mag = np.hypot(dx, dy)
        good = mag > 1e-12
        x0, y0, x1, y1 = x0[good], y0[good], x1[good], y1[good]
        dx, dy, mag = dx[good], dy[good], mag[good]

        if len(x0) == 0:
            self.clear()
            return

        ux, uy = dx / mag, dy / mag
        head = 0.28 * arrow_len
        nan = np.full(len(x0), np.nan)

        sx = np.stack([x0, x1, nan], axis=1).ravel()
        sy = np.stack([y0, y1, nan], axis=1).ravel()

        lx = x1 - head * (ux + uy)
        ly = y1 - head * (uy - ux)
        rx = x1 - head * (ux - uy)
        ry = y1 - head * (uy + ux)

        hx = np.stack([x1, lx, nan, x1, rx, nan], axis=1).ravel()
        hy = np.stack([y1, ly, nan, y1, ry, nan], axis=1).ravel()

        self._shaft.setData(sx, sy)
        self._head.setData(hx, hy)
        