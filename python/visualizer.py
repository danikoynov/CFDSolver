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
    
    def __init__(self, sim, grid_w, grid_h, cell_size):
        self.sim = sim
        self.grid_w = grid_w
        self.grid_h = grid_h
        self.cell_size = cell_size
        
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
    def draw_arrow(screen, color, start, vec, max_len=16):
        vx, vy = vec
        mag = (vx * vx + vy * vy) ** 0.5
        if mag < 1e-10:
            return

        scale = min(max_len / mag, 1e9)
        dx = vx * scale
        dy = vy * scale

        x0, y0 = start
        x1 = x0 + dx
        y1 = y0 + dy

        pygame.draw.line(screen, color, (x0, y0), (x1, y1), 2)

        length = (dx * dx + dy * dy) ** 0.5
        ux = dx / length
        uy = dy / length
        ah = max(4, 0.25 * length)

        left = (x1 - ah * (ux + uy), y1 - ah * (uy - ux))
        right = (x1 - ah * (ux - uy), y1 - ah * (uy + ux))

        pygame.draw.line(screen, color, (x1, y1), left, 2)
        pygame.draw.line(screen, color, (x1, y1), right, 2)
        
    def draw_pressure_and_velocity(self, screen, grid, cell_size, font):
        pressures = self.collect_pressure(grid)
        vf = grid.velocity()
        w = grid.width()
        h = grid.height()

        flat = [p for row in pressures for p in row]
        pmin = min(flat)
        pmax = max(flat)

        screen.fill((20, 20, 20))

        for y in range(h):
            for x in range(w):
                p = pressures[y][x]
                color = self.pressure_to_color(p, pmin, pmax)

                screen_y = (h - 1 - y) * cell_size
                rect = pygame.Rect(x * cell_size, screen_y, cell_size, cell_size)
                pygame.draw.rect(screen, color, rect)
                pygame.draw.rect(screen, (40, 40, 40), rect, 1)

        # draw fewer arrows on dense grids
        stride = 1
        if max(w, h) <= 20:
            stride = 1
        elif max(w, h) <= 35:
            stride = 2
        elif max(w, h) <= 60:
            stride = 3
        else:
            stride = 4

        max_speed = 1e-12
        for y in range(0, h, stride):
            for x in range(0, w, stride):
                u, v = self.cell_center_velocity(vf, x, y)
                speed = (u * u + v * v) ** 0.5
                max_speed = max(max_speed, speed)

        arrow_spacing = stride * cell_size
        arrow_len = max(8, 0.45 * arrow_spacing)

        for y in range(0, h, stride):
            for x in range(0, w, stride):
                u, v = self.cell_center_velocity(vf, x, y)

                sx = (x + 0.5) * cell_size
                sy = (h - y - 0.5) * cell_size

                vx = (u / max_speed) * arrow_len
                vy = -(v / max_speed) * arrow_len

                self.draw_arrow(screen, (255, 255, 255), (sx, sy), (vx, vy), max_len=arrow_len)

        text = font.render(
            f"{w}x{h} | Pressure min={pmin:.4f} max={pmax:.4f}",
            True,
            (255, 255, 255),
        )
        screen.blit(text, (10, 10))

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
    def run(self):
        window_width = self.grid_w * self.cell_size
        window_height = self.grid_h * self.cell_size

        pygame.init()
        screen = pygame.display.set_mode((window_width, window_height))
        pygame.display.set_caption("CFD Pressure + Velocity Visualization (30x30)")
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

            if not paused:
                self.sim.tick()

            self.print_velocity_field()
            self.draw_pressure_and_velocity(screen, self.sim.grid(), self.cell_size, font)
            clock.tick(30)

        pygame.quit()        
        
    