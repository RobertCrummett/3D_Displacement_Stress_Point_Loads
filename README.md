# Displacements and Stress Due to Point Loads
----------------------------------------------
This file contains functions to calculate the displacements and stress in the ground
as a result of a point force at the origin.

----------------------------------------------
Functions Displacement_3D_NormalPointLoad() and Displacement_3D_TangentialPointLoad()
solve for displacements induced by a normal and tangential (+x direction) force 
respectively.

Functions Traction_3D_NormalPointLoad() and Traction_3D_TangentialPointLoad() solve
for tractions (stress) induced by a normal and tangential (+x direction) force
respectively.

Function Interpolate_Traction_3D takes outputs of Traction_3D_NormalPointLoad() or
Traction_3D_TangentialPointLoad() and interpolated the meshes, giving final plots and
graphics a smoother look.

The script Examples gives an example of how to use the functions and define their inputs.

---------------------------------------------
I have undated this repository with the code Traction_3D_NormalLoad(). For a load shape
that can be characerized in terms of [x_min,x_max] and [y_min,y_max], this code will
integrate the point load solutions to produce a surface traction.

# Sources
----------------------------------------------
The source for this work is:

Jaeger, J.C., Cook, N.G.W. and Zimmerman, R.W. (2007) Fundamentals of Rock 
Mechanics. 4th Edition, Chapman and Hall, London. pg. 408-411

Thank you to Dr. Hammond for loaning me the book.
