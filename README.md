
Code simulates two-body, three-body and N-body gravity. 

Please run the code cell-by-cell, otherwise the animations will all run at once and cause code to
slow dramatically. Best way to run in VSCode is using the jupyter notebook file (.ipynb) as 
notebooks has better implemetation of running cell-by-cell. Can run natively in VSCode python but 
not as fast. Animations run smoothest in Spyder.

IMPORTANT: IF RUNNING IN SPYDER OR JUPYTER NOTEBOOK (INCLUDING NOTEBOOK INSIDE VSCODE), THEN 
THE LINE '%matplotlib qt' MUST BE INCCLUDED WITH THE LIBRARY IMPORTS. 
    THIS ENABLES INTERACTIVE FIGURES AND ANIMATIONS. 

Cell 1:
Imports all modules/libraries needed. 

Cell 2: 
Sets inital conditions and solves a two-body 2D system of Earth and the Sun.

Cell 3:
Plots distance of Earth from Sun over time. (Test figure 1).

Cell 4:
Plots orbit view from above. X and Y solutions are plotted for Earth and the Sun. 

Cell 5: 
Repeats the plot in cell 4 but in a 3D axis. Sets number of ticks on axis to 6 to avoid
overlapping numbers.

Cell 6: 
Sets inital conditions of Mars, then solves the system of Sun and Mars. 
Then adds the two-body solution of Mars to a 3D plot of Earth and the Sun.

Cell 7: 
Solves three bodies in 3D. In this case we use two heavy stars and one smaller mass (planet). 
Then animates a line plot of position for all objects. Use RK23 for circular orbits and RK45 for 
orbits that are more chaotic. 

Cell 8: 
Plots x,y,x coordinates of the previous animation over time. 
Repeats the plot in cell 4 but in a 3D axis. Sets number of ticks on axis to 6 to avoid overlapping numbers.
Pad used in both label and ticks to ensure axes still readable. 

Cell 9: 
Solves another 3D three body problem, this time one star and two planets. Sets x and y positions of all objects
to a line object then animates the lines. 

Cell 10: 
Sets up functions to incrimentally calculate force and acceleration. Then a function to animate the N-body
line. 

    Cell 10 Defines:

    -n = the number of bodies,

    -seed number,

    -vx,vy,vz = random velocities in each direction using random uniform distribution. 

    -positions of objects in each direction using random uniform distribution. 

    -masses= random masses of each object uniformly distributed between two limits in terms of solar masses. 


