import os
import sys
from pathlib import Path

root = Path(__file__).resolve().parent.parent
build_dir = root / "build"

if hasattr(os, "add_dll_directory"):
    os.add_dll_directory(r"C:\msys64\ucrt64\bin")

sys.path.insert(0, str(build_dir))

import cfdsolver_py
import pygame


class Visualizer:

    SOLID_COLOR = (120, 120, 120)
    BACKGROUND_COLOR = (20, 20, 20)
    ARROW_COLOR = (255, 255, 255)
    TEXT_COLOR = (0, 200, 230)
    GRID_LINE_COLOR = (80, 80, 80)  # thin grid line color

    def __init__(
        self,
        sim,
        grid_w,
        grid_h,
        cell_size,
        scale_arrows_by_magnitude=True,
    ):
        self.sim = sim
        self.grid_w = grid_w
        self.grid_h = grid_h
        self.cell_size = cell_size
        self.scale_arrows_by_magnitude = scale_arrows_by_magnitude

    def collect_pressure(self, grid):
        pf = grid.pressure()
        w = grid.width()
        h = grid.height()

        values = []
        for y in range(h):
            row = []
            for x in range(w):
                row.append(pf.get(x, y))
            values.append(row)

        return values

    def collect_fluid_pressures(self, grid, pressures):
        bc = grid.boundary_conditions()
        w = grid.width()
        h = grid.height()

        values = []
        for y in range(h):
            for x in range(w):
                if bc.type(x, y) != cfdsolver_py.CellType.SOLID:
                    values.append(pressures[y][x])

        return values

    @staticmethod
    def pressure_to_color(p, pmin, pmax):
        if abs(pmax - pmin) < 1e-12:
            t = 0.5
        else:
            t = (p - pmin) / (pmax - pmin)

        t = max(0.0, min(1.0, t))

        r = int(255 * t)
        g = 0
        b = int(255 * (1.0 - t))

        return (r, g, b)

    @staticmethod
    def cell_center_velocity(vf, x, y):
        u = 0.5 * (vf.get_u(x, y) + vf.get_u(x + 1, y))
        v = 0.5 * (vf.get_v(x, y) + vf.get_v(x, y + 1))
        return u, v

    @staticmethod
    def draw_arrow(screen, color, start, vec, max_len=None):
        vx, vy = vec
        mag = (vx * vx + vy * vy) ** 0.5

        if mag < 1e-10:
            return

        if max_len is not None and mag > max_len:
            scale = max_len / mag
            vx *= scale
            vy *= scale
            mag = max_len

        x0, y0 = start
        x1 = x0 + vx
        y1 = y0 + vy

        pygame.draw.line(screen, color, (x0, y0), (x1, y1), 2)

        ux = vx / mag
        uy = vy / mag
        arrow_head = max(4, 0.25 * mag)

        left = (x1 - arrow_head * (ux + uy), y1 - arrow_head * (uy - ux))
        right = (x1 - arrow_head * (ux - uy), y1 - arrow_head * (uy + ux))

        pygame.draw.line(screen, color, (x1, y1), left, 2)
        pygame.draw.line(screen, color, (x1, y1), right, 2)

    def draw_grid_lines(self, screen, w, h, cell_size):
        # vertical lines
        for x in range(w + 1):
            px = x * cell_size
            pygame.draw.line(
                screen,
                self.GRID_LINE_COLOR,
                (px, 0),
                (px, h * cell_size),
                1,
            )

        # horizontal lines
        for y in range(h + 1):
            py = y * cell_size
            pygame.draw.line(
                screen,
                self.GRID_LINE_COLOR,
                (0, py),
                (w * cell_size, py),
                1,
            )

    def draw_pressure_and_velocity(self, screen, grid, cell_size, font):
        pressures = self.collect_pressure(grid)
        vf = grid.velocity()
        bc = grid.boundary_conditions()
        w = grid.width()
        h = grid.height()

        fluid_pressures = self.collect_fluid_pressures(grid, pressures)
        if fluid_pressures:
            pmin = min(fluid_pressures)
            pmax = max(fluid_pressures)
        else:
            pmin = 0.0
            pmax = 0.0

        screen.fill(self.BACKGROUND_COLOR)

        for y in range(h):
            for x in range(w):
                screen_y = (h - 1 - y) * cell_size
                rect = pygame.Rect(x * cell_size, screen_y, cell_size, cell_size)

                if bc.type(x, y) == cfdsolver_py.CellType.SOLID:
                    color = self.SOLID_COLOR
                else:
                    p = pressures[y][x]
                    color = self.pressure_to_color(p, pmin, pmax)

                pygame.draw.rect(screen, color, rect)

        # draw thin grid lines on top of cells
        self.draw_grid_lines(screen, w, h, cell_size)

        stride = 1
        if max(w, h) <= 20:
            stride = 1
        elif max(w, h) <= 35:
            stride = 2
        elif max(w, h) <= 60:
            stride = 3
        else:
            stride = 4

        arrow_spacing = stride * cell_size
        arrow_len = max(8, 0.45 * arrow_spacing)

        max_speed = 1e-12
        for y in range(0, h, stride):
            for x in range(0, w, stride):
                if bc.type(x, y) == cfdsolver_py.CellType.SOLID:
                    continue

                u, v = self.cell_center_velocity(vf, x, y)
                speed = (u * u + v * v) ** 0.5
                max_speed = max(max_speed, speed)

        for y in range(0, h, stride):
            for x in range(0, w, stride):
                if bc.type(x, y) == cfdsolver_py.CellType.SOLID:
                    continue

                u, v = self.cell_center_velocity(vf, x, y)
                speed = (u * u + v * v) ** 0.5

                if speed < 1e-12:
                    continue

                sx = (x + 0.5) * cell_size
                sy = (h - y - 0.5) * cell_size

                if self.scale_arrows_by_magnitude:
                    scale = arrow_len / max_speed
                    vx = u * scale
                    vy = -v * scale
                else:
                    vx = arrow_len * u / speed
                    vy = -arrow_len * v / speed

                self.draw_arrow(
                    screen,
                    self.ARROW_COLOR,
                    (sx, sy),
                    (vx, vy),
                    max_len=arrow_len,
                )

        text = font.render(
            f"{w}x{h} | Pressure min={pmin:.4f} max={pmax:.4f}",
            True,
            self.TEXT_COLOR,
        )
        screen.blit(text, (10, 10))

        mode_text = font.render(
            f"Arrow mode: {'scaled' if self.scale_arrows_by_magnitude else 'fixed'}",
            True,
            self.TEXT_COLOR,
        )
        screen.blit(mode_text, (10, 35))

        pygame.display.flip()

    def print_velocity_field(self):
        vf = self.sim.grid().velocity()
        w = self.grid_w
        h = self.grid_h

        print("\nVelocity field (cell-centered):")
        for y in range(h - 1, -1, -1):
            row = []
            for x in range(w):
                u, v = self.cell_center_velocity(vf, x, y)
                row.append(f"({u: .3f},{v: .3f})")
            print(" ".join(row))

    def print_pressure_field(self):
        pf = self.sim.grid().pressure()
        w = self.grid_w
        h = self.grid_h

        print("\nPressure field:")
        for y in range(h - 1, -1, -1):
            row = []
            for x in range(w):
                row.append(f"{pf.get(x, y):7.3f}")
            print(" ".join(row))

    def run(self):
        window_width = self.grid_w * self.cell_size
        window_height = self.grid_h * self.cell_size

        pygame.init()
        screen = pygame.display.set_mode((window_width, window_height))
        pygame.display.set_caption("CFD Pressure + Velocity Visualization")
        font = pygame.font.SysFont(None, 24)
        clock = pygame.time.Clock()

        running = True
        paused = False

        while running:
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    running = False
                elif event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_SPACE:
                        paused = not paused
                    elif event.key == pygame.K_a:
                        self.scale_arrows_by_magnitude = (
                            not self.scale_arrows_by_magnitude
                        )

            if not paused:
                self.sim.tick()

            self.print_pressure_field()
            self.print_velocity_field()
            self.draw_pressure_and_velocity(
                screen,
                self.sim.grid(),
                self.cell_size,
                font,
            )
            clock.tick(30)

        pygame.quit()