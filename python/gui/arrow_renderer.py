import numpy as np
import pyqtgraph as pg


class ArrowRenderer:
    def __init__(self, plot: pg.PlotWidget):
        self.plot = plot
        self.items: list[pg.PlotDataItem] = []

    def clear(self) -> None:
        for item in self.items:
            self.plot.removeItem(item)

        self.items.clear()

    def draw(
        self,
        u: np.ndarray,
        v: np.ndarray,
        solid: np.ndarray,
        density: int,
        uniform_length: bool,
    ) -> None:
        self.clear()

        h, w = u.shape
        side = max(w, h)

        density = max(1, density)
        stride = max(1, side // density)

        speed_all = np.hypot(u, v)
        max_speed = float(speed_all.max())

        if max_speed < 1e-12:
            return

        js = np.arange(0, h, stride)
        is_ = np.arange(0, w, stride)

        jj, ii = np.meshgrid(js, is_, indexing="ij")
        mask = ~solid[jj, ii]

        jj = jj[mask]
        ii = ii[mask]

        ui = u[jj, ii]
        vi = v[jj, ii]

        speed = np.hypot(ui, vi)
        valid = speed > 1e-12

        jj = jj[valid]
        ii = ii[valid]
        ui = ui[valid]
        vi = vi[valid]
        speed = speed[valid]

        if len(ii) == 0:
            return

        x0 = ii + 0.5
        y0 = jj + 0.5

        arrow_len = 0.42 * stride

        if uniform_length:
            dx = ui / speed * arrow_len
            dy = vi / speed * arrow_len
        else:
            dx = ui / max_speed * arrow_len
            dy = vi / max_speed * arrow_len

        x1 = x0 + dx
        y1 = y0 + dy

        mag = np.hypot(dx, dy)
        good = mag > 1e-12

        x0 = x0[good]
        y0 = y0[good]
        x1 = x1[good]
        y1 = y1[good]
        dx = dx[good]
        dy = dy[good]
        mag = mag[good]

        if len(x0) == 0:
            return

        ux = dx / mag
        uy = dy / mag

        head = 0.28 * arrow_len

        lx = x1 - head * (ux + uy)
        ly = y1 - head * (uy - ux)

        rx = x1 - head * (ux - uy)
        ry = y1 - head * (uy + ux)

        nan = np.full(len(x0), np.nan)

        sx = np.stack([x0, x1, nan], axis=1).ravel()
        sy = np.stack([y0, y1, nan], axis=1).ravel()

        hx = np.stack([x1, lx, nan, x1, rx, nan], axis=1).ravel()
        hy = np.stack([y1, ly, nan, y1, ry, nan], axis=1).ravel()

        pen = pg.mkPen((180, 190, 220, 200), width=1)
        pen_head = pg.mkPen((180, 190, 220, 140), width=1)

        shaft_item = pg.PlotDataItem(sx, sy, pen=pen, connect="finite")
        head_item = pg.PlotDataItem(hx, hy, pen=pen_head, connect="finite")

        self.plot.addItem(shaft_item)
        self.plot.addItem(head_item)

        self.items.extend([shaft_item, head_item])
        