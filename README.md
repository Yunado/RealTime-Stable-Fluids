# Physics-Based Animation: 2D Real Time Stable Fluids
2D Stable Fluids
Python implementation of Stable Fluids by Jos Stam with Numpy and Pygame.

## Execution
Makesure to have Numpy and Pygame installed in ur Python environment.
Note that on the first time installing Pygame, you may require to restart the computer,
or else MEM allocation error from Pygame will prompt.

Makesure to put main.py, simulation.py and stable_fluids.py inside one folder.
Then start the program by running main.py.
parameters are at the top of main.py, from line 6 to 10.
Diss = Dissipation, Diff = Diffusion, VISC = Viscosity, DT = dt, N = size of the grid.
The default parameters are set to simulate Gaseous phenomena similar to the Stable Fluids paper.

## Interaction
After start the GUI, press down the left mouse button to put in sourses, 
where drag the mouse will create velocity.
Press "S" button to start the simulation once you draw on the canvas.
Press "A" button to pause the simulation.
Press "R" button will reset the simulation.

## Video
https://www.youtube.com/watch?v=wprjBG5b5VA
