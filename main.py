from simulation import *
import pygame
import math
import random

DISS = 0.01
DIFF = 0.0001
VISC = 0.0001
DT = 0.1
N = 46
SIZE = N + 2
SCALE = 10
FPS = 24

def setup():
    simulation = Simulation(DT, N, DISS, DIFF, VISC)
    return simulation

def main():
    sim = setup()
    pygame.init()
    screen = pygame.display.set_mode((SIZE * SCALE, SIZE * SCALE))
    clock = pygame.time.Clock()
    pygame.display.set_caption("STABLE FLUID")
    running = True
    hold = False
    sim_start = False
    pos = pygame.mouse.get_pos()

    def render_density():
        for i in range(SIZE):
            for j in range(SIZE):
                d = sim.fluid.d[sim.fluid.ind(i, j)]
                if d > 255:
                    d = 255
                if d < 0:
                    d = 0
                pygame.draw.rect(screen, [d, d, d], [j * SCALE, i * SCALE, SCALE, SCALE])

    # def render_velocity():
    #     for i in range(SIZE):
    #         for j in range(SIZE):
    #             d = sim.fluid.d[sim.fluid.ind(i, j)]
    #             if d > 255:
    #                 d = 255
    #             if d < 0:
    #                 d = 0
    #             u = math.floor(sim.fluid.Vx[sim.fluid.ind(i, j)])
    #             v = math.floor(sim.fluid.Vy[sim.fluid.ind(i, j)])
    #             theta = 0
    #             if v > 0 and u > 0:
    #                 theta = min(math.degrees(math.atan(u/v)), math.degrees(math.atan(v/u)))
    #             uv = pygame.math.Vector2(u, v)
    #             start_p = pygame.math.Vector2(j * SCALE, i * SCALE)
    #             end_p = start_p
    #             if uv.length() > 0:
    #                 end_p = start_p + uv.rotate(theta)
    #             # enb_p = start_p + end_p.rotate(theta)
    #             pygame.draw.line(screen, [d/2, d, d], start_p, end_p)

    while running:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
            if event.type == pygame.MOUSEBUTTONDOWN and event.button == 1:
                hold = True
                pos = pygame.mouse.get_pos()
            if event.type == pygame.MOUSEBUTTONUP and event.button == 1:
                hold = False
            if hold:
                curr_pose = pygame.mouse.get_pos()
                diffx = curr_pose[0] - pos[0]
                diffy = curr_pose[1] - pos[1]
                print(curr_pose)
                print(diffx, diffy)
                sim.fluid.add_density(math.floor(curr_pose[1]/SCALE), math.floor(curr_pose[0]/SCALE), random.randint(100, 200))
                sim.fluid.add_velocity(math.floor(curr_pose[1]/SCALE), math.floor(curr_pose[0]/SCALE), diffx, diffy)
            if event.type == pygame.KEYDOWN:
                # press s to start the simulation
                if event.key == pygame.K_s:
                    sim_start = True
                # press r to reset the simulation
                if event.key == pygame.K_r:
                    sim.reset()
                    sim_start = False
                # press a to pause the simulation
                if event.key == pygame.K_a:
                    sim_start = False

        if sim_start:
            sim.simulate()
        render_density()
        # density_2d = sim.fluid.d.reshape(SIZE, SIZE)
        # surface = pygame.surfarray.make_surface(density_2d)
        # surface = pygame.transform.scale(surface, (SIZE * SCALE, SIZE * SCALE))
        # screen.blit(surface, (0, 0))
        pygame.display.update()
        clock.tick(FPS)


if __name__ == "__main__":
    main()
    pygame.quit()
