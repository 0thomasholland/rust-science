# Critical Cellular Automation

Written in Fortran, this code is based on a 3D grid that has the following conditions:

- Each cell (i,j,k) on the grid has a positive integer value
- Each cell on the grid has a critical maximum of a value of 6; above this value it redistributes it to the surrounding 6 cells (e.g. i+1,j+1,k+1 += 1) and sets it value to what it was less 6 (i,j,k -= 6).

Initial starting conditions:

- Either from a blank canvas
- Random values between 0-5 randomly selected and allocated to each grid point

Within each logic loop:

- Checks if any of the cells are above 6
- If they are above 6 the points redistribute
- This check and redistribute occurs until reaches an equilibrium, and each time through the state is written to disk via appending a file
- Then a random point is chosen and 1 is added to it's value

A python script accompanies which takes the output file and plots a 3D revolving video of this, where the opacity of each 3D cell is related to its criticality

## Python Visualization

Needs to have numpy and matplotlib pip installed and ffmpeg also installed.

```bash
python3 visualize.py --output animation.mp4 --fps 15 --grid-size 20 20 20 
```
