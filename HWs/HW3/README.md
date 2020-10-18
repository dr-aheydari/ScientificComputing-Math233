# HW3: Simulating Simple Advection in 2D with Semi-Lagrangian and LevelSets
This repository will host my code for solving 2D advection equation using Semi-Lagrangian method, and also using Level Sets with reinitialization. This repository is meant to accompany my LaTeX write-up turned in online. Please be sure to ask permission/cite this work if you use this for any graded assignment. 

## Semi-Lagrangian: Unconditional Stability
 
In this section, we show simulations of the SL method on a wide range of time steppings that often cause instability due to the CFL condition. For faster computation time, all visualizations are shown on a 128X128 grid, although this code has been tested for uniform grids up to 512X512.

 
### dt = 0.25(dx)
The simulation below has been sped up to 8x the actual evolution time. For the actual evolution simulation, please refer to the linked Box folder below: 

![dt=0.25dx](https://media.giphy.com/media/G1sAoEurj3sf5SrxHY/giphy.gif)

[dt=0.25dx Simulation Video (Downloadable)](https://ucmerced.box.com/s/568fr7rcbmxjx9syqr45k1dz2sa71rzt)


### dt = 0.5 dx
The visualization below is 6x the speed of the actual evolution time.


![dt=0.5dx GIF](https://media.giphy.com/media/nVzPMWdtvl1AdOethD/giphy.gif)

[dt=0.5dx Simulation Video (Downloadable)](https://ucmerced.box.com/s/568fr7rcbmxjx9syqr45k1dz2sa71rzt)

### dt = dx 
The visualization below is 4x the speed of the actual evolution time. Original evolution time can be found in the linked Box folder below.

![dt = dx GIF](https://media.giphy.com/media/508XKw7hdg8PoZBina/giphy.gif)


[dt = dx Video (Downloadable)](https://ucmerced.box.com/s/568fr7rcbmxjx9syqr45k1dz2sa71rzt)


### dt = 2dx
The visualization illustrates the actual evolution time.

![dt = 2dx GIF](https://media.giphy.com/media/508XKw7hdg8PoZBina/giphy.gif)


[dt = 2dx Video (Downloadable)](https://ucmerced.box.com/s/568fr7rcbmxjx9syqr45k1dz2sa71rzt)

### dt = 5dx

The visualization illustrates the actual evolution time. 

![dt = 5dx GIF](https://media.giphy.com/media/508XKw7hdg8PoZBina/giphy.gif)

[dt = 5dx Video (Downloadable)](https://ucmerced.box.com/s/568fr7rcbmxjx9syqr45k1dz2sa71rzt)

### dt = 10dx
The visualization illustrates the actual evolution time. 

![dt = 10dx GIF](https://media.giphy.com/media/508XKw7hdg8PoZBina/giphy.gif)

[dt = 10dx Video (Downloadable)](https://ucmerced.box.com/s/568fr7rcbmxjx9syqr45k1dz2sa71rzt)


## Level Set: Testing Reinitialization 
