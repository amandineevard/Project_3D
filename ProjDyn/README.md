# Projective Dynamics and ShapeUp Framework

## What it is

A minimal framework that implements Projective Dynamics simulations or ShapeUp style mesh editing.
The GUI is implemented using **nanogui**, PD and ShapeUp are implemented using **Eigen** and **OpenMP**.
The mesh is intermediately loaded as a **Surface_mesh** object, but converted to simple vertex positions / triangle indices matrices for PD and ShapeUp internally.
**tetgen** is used to generate tetrahedral meshes if volumetric constraints are desired.

The available constraints are: position constraints, edge springs, triangle strain and bending, tetrahedral strain, floor collisions.

## Build
### Unix
This project uses Wenzel Jakob's [nanogui](https://github.com/wjakob/nanogui), which is included under externals,
but has additional dependencies, which can be resolved on Ubuntu/Debian, via:

    $ apt-get install cmake xorg-dev libglu1-mesa-dev

and on Redhead/Fedora, via:

	$ sudo dnf install cmake mesa-libGLU-devel libXi-devel libXcursor-devel libXinerama-devel libXrandr-devel xorg-x11-server-devel

After that, simply build the project as usual:

	mkdir build
	cd build
	cmake ..
	make

    
To run the ShapeUp editor, enter the build directory and run

    ./shapeup/shapeup
    
To run the Projective Dynamics simulation, enter the build directory and run

	./simulation/simulation

### Windows/MVSC 19
	Run CMake
	Set the source code dir to the base dir
	Set the build dir to the base dir + "/build"
	Press Configure
	Press Generate
	Open projdyn.sln in the build directory
	Build solution and run the simulation project or the shapeup project

## Simulation
	Open a mesh.
	Press Tetrahedralize if tets are desired (e.g. Tet strain constraints).
	Add Constraints (e.g. Floor constraints and Tet strain).
	Press Run Simulation.
	Constraints can be added/removed in a running simulation.
	The camera can be controlled via left-click + drag to rotate, or right click to translate.
	The mesh can be grabbed and moved via alt/option + left-click + drag.
	If anything goes wrong, use Reset Positions and Clear Constraints and try again.
    
## ShapeUp
    Open a mesh.
    Vertices can be selected using ctrl + left-click + drag, which shows a selection box.
    Selected vertices (blue) can be made into position constraint groups using the button “Fix Selection”. Create two such groups on different places on the mesh.
    Vertices marked in green can be moved (as a group) by alt/option + left-click + drag, which changes the target positions of the constraints.
    The shape is maintained by using either edge springs or triangle strain (or, once implemented, Iso 1-Ring constraints), which can be selected using the corresponding radio button. Note how much nicer the shape is maintained using the triangle strain + bending constraints, see the next slides for details.
    Additional constraints can be added to selected vertices by the buttons to the lower left. For now, only smoothing constraints are available.
    The currently used constraints are shown as groups in the window to the top right. There, the constraint weights can be changed using the sliders. Hovering over a slider, shows the vertices on which these constraints are defined.