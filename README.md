
<h1 style="color:teal;">Finite volume solver for shallow water moment solver</h1>


This Github repository contains the source files of a solver written in Fortran90 tailored to approximate a variety of shallow water moment models in one or two spatial dimensions. 
The development of this software is part of the project HiWAVE - Natural hazard prediction with adaptive hierarchical wave models, of the research programme Vidi ENW (https://doi.org/10.61686/CBVAB59929).

The benchmark partial differential equation (PDE), in the one-dimensional case,  solved by this software is the following first-order system of $`n`$ equations:
```math
\color{brown}\dfrac{\partial\boldsymbol{U}}{\partial t} \,+\, \dfrac{\partial \boldsymbol{F}}{\partial x}(\boldsymbol{U}) \,+\, \boldsymbol{A}(\boldsymbol{U})\dfrac{\partial \boldsymbol{U}}{\partial x} \,=\, \boldsymbol{B}(\boldsymbol{U})\dfrac{\partial \boldsymbol{U}}{\partial x} \,+\, \boldsymbol{S}(\boldsymbol{U}),\qquad \qquad(1)
```
where $`t>0`$ is time, $`x`$ is the spatial coordinate, $`\boldsymbol{U}\in \mathbb{R}^{n}`$ is the vector of unknowns, $`\boldsymbol{F}`$ is a flux, $`\boldsymbol{A}`$ is a nonlinear matrix function which is typically the Jacobian of a vector function, $`\boldsymbol{B}`$ is a remaining non-conservative matrix, and  $`\boldsymbol{S}`$ is the source term. An alternative form of writting system (1) is the non-conservative compact form:
```math
\color{brown}\dfrac{\partial\boldsymbol{U}}{\partial t}  \,+\, \Big( \partial_{\boldsymbol{U}}\boldsymbol{F}(\boldsymbol{U}) + \boldsymbol{A}(\boldsymbol{U}) - \boldsymbol{B}(\boldsymbol{U}) \Big) \dfrac{\partial \boldsymbol{U}}{\partial x} \,=\,  \boldsymbol{S}(\boldsymbol{U}).\qquad \qquad(2)
```
The numerical scheme employed to solve the PDE combines a variety of ingredients:
- High-order time approximations for $`\partial_t \boldsymbol{U}`$.
- Conservative flux approximation related to the conservative part $`\boldsymbol{F}`$ denoted by $`F_{j+1/2}`$.
- Generalized Roe matrices for the non-conservative terms $`\boldsymbol{A}`$ and $`\boldsymbol{B}`$, or for the combinations $`\boldsymbol{A}-\boldsymbol{B}`$ or $`\partial_{\boldsymbol{U}}\boldsymbol{F}+ \boldsymbol{A} - \boldsymbol{B}`$. The generalized Roe matrices are computed via path conservative integrals for the fluctuations denoted by $`D_{j+1/2}^+`$ and  $`D_{j+1/2}^-`$.
- Explicit or implicit treatment of the source term via splitting and using Newton-Raphson solver.
- High-order polynomial reconstructions and slope limiters.

The components of the vector of unknown $`\boldsymbol{U}`$ are called the conservative variables. A second vector of unknown used throughout the subroutines of the program is the vector of primitive components $`\boldsymbol{V}`$ which is defined as follows
```math
V_1 = U_1 \quad\text{and}\quad V_i = U_i/U_1,\qquad \forall i=2,\dots,n
```

Locally at each subroutine, the conservative vector $`\boldsymbol{U}`$ is denoted as Cons and the primitive vector $`\boldsymbol{V}`$ is denoted as Prim. The conversion between the primitive variables to the conservative and viceversa is made trhough the subroutines PrimToCons and ConsToPrim, respectively.

This program is primarily designed for **Linux** operative system. A **Windows** version can be done by modifying the make file in the main folder and adapting the lines of code related to generation of files and folders management. 

### Repositories and organization

The software is organized in three subfolders:
- **src**: Contains all the f90 source codes, including the additional files **deplist.mk**, **sources.mk** and **makefile**, which are required for the compilation of the main program.
- **lib**: This is a repository where compressed auxiliary files arising only after compilation will be placed. A repository not used for further purposes but needs to be there for the program to compile.
- **bin**: Folder where the executable program called **main** will be created. This is the folder where the examples (located in subfolders) will be run.

The file called **makefile** at the main shallow-water-moments repository (different than the one in src) is the one used to compile the code in the terminal.

## Compilation
To compile the program, the default option is to use the GNU compiler **gfortran** in **Linux**, which can be installed running the following lines in the terminal:

```console
sudo apt update
sudo apt install gfortran
```

Then at the main folder *shallow-water-moments* simply run the make file by typing in the terminal:
```console
user@pc:~/shallow-water-moments$ make
```
The executable program **main** will be created at the folder *bin*, several ".o" and ".mod" files are going to be created in the folder *src* and an additional file called libmain.a is also generated in the folder *lib*. To clean the compilation, the following line can be run:
```console
user@pc:~/shallow-water-moments$ make clean
```
Then all files generated after compilation are deleted.

### Makefile:

There are two makefile, one in the main folder *shallow-water-moments* that we denote here by makefile 1, which is the one to be executed in the terminal, and one in the subfolder *src* denoted by makefile 2, which is called by the first makefile. The makefile 2 uses the auxiliary files deplist.mk (structure of the modules used by each f90 file) and sources.mk (list of f90 files included in the program). In the fist part of makefile 2, the compiler is chosen between gfortran and intel, and flags to set: 

- language specification
- float number precision
- backslash statement
- number of characters per line
- setting floating point exception
- settings for debugging

The makefile also contains an additional flag for using **Blas** and **Lapack** (Fortran libraries) which are required to compile the program:

```console
BLAS_LAPACK = -llapack -lblas
```

To run the two dimensional examples available in the software the following line needs to be uncommented:

```console
#FCFLAGS += -I. -DTWODIM
```

## Code structure
The f90 containing the program and all modules are in the source folder *src*. The file **main.f90** comprises the program **FiniteVolumeSWM** and therefore which contains the main structure of the code including the initialization, mesh generation, time loop, saving and postprocessing procedures. The structure of the main file is displayed in the following diagram:

<img width="1200" height="900" alt="DiagramFVSWME" src="https://github.com/user-attachments/assets/2559f9c5-d8db-47ce-b151-b525093a0220" />


The modules and subroutines called by this program are:

- ***finitevolume_vars.f90*** (module $\color{teal}\texttt{MOD\\_FiniteVolume}$): In this module the majority of the global variables used throughout the many subroutines are defined. It also includes variables that are allocated in further subroutines. This module is used in almost all subroutines of the program. 
- ***parameters.f90*** (module $\color{teal}\texttt{MOD\\_Parameters}$): The purpose of this module is to read the input file, define parameters and set the pointers depending on the flags given at the input file.
 
   $\quad\bullet$ $\color{blue}\texttt{Current 1D models:}$<br>
   $\qquad\diamond$ Shallow water equations (SWE)<br>
   $\qquad\diamond$ Shallow water moment equations (SWME)<br>
   $\qquad\diamond$ Shallow water linearized moment equations (SWLME)<br>
   $\qquad\diamond$ Hyperbolic SWME (HSWME)<br>
   $\qquad\diamond$ Inclined HSWME<br>
   $\qquad\diamond$ Beta hyperbolic SWME (HSWME)<br>
   $\qquad\diamond$ Moment regularization for hyperbolic SWME (MHSWME)<br>
   $\qquad\diamond$ Primitive regularization for hyperbolic SWME (PHSWME)<br>
   $\qquad\diamond$ Primitive moment regularization for hyperbolic SWME (PMHSWME)<br>

   $\quad\bullet$ $\color{blue}\texttt{Friction models:}$<br>
   $\qquad\diamond$ Newtonian slip<br>
   $\qquad\diamond$ Newtonian Manning<br>
   $\qquad\diamond$ Coulomb type<br>
   $\qquad\diamond$ Granular<br>
   $\qquad\diamond$ Savage Hutter<br>
   
   $\quad\bullet$ $\color{blue}\texttt{FVM schemes:}$<br>
   $\qquad\diamond$ Path intergals with: Linear, Quadratic or Power law paths<br>
   $\qquad\diamond$ Viscosity matrices: Lax-Friedrichs, Lax Wendroff, Force (Price-C), HLL<br>
   
   $\quad\bullet$ $\color{blue}\texttt{Time integrators:}$<br>
   $\qquad\diamond$ Forward Euler<br>
   $\qquad\diamond$ SSPRK4<br>
   $\qquad\diamond$ RK65<br>
   $\qquad\diamond$ Alternative implicit solver for the source term (Nonlinear solver)<br>

   $\quad\bullet$ $\color{blue}\texttt{Time integrators:}$<br>
   $\qquad\diamond$ Periodic<br>
   $\qquad\diamond$ Transmissive<br>
   $\qquad\diamond$ Inflow<br>
   $\qquad\diamond$ Outflow<br>
   $\qquad\diamond$ Reflecting<br>

   $\quad\bullet$ $\color{blue}\texttt{Partially available (under development)}$<br>
   $\qquad\diamond$ Fully 2D problems<br>
   $\qquad\diamond$ Eigen value method (upwind) for the approximation of the linearized Roe matrix<br>
   $\qquad\diamond$ Riemann solver for conservative fluxes<br>
   $\qquad\diamond$ Polynomial reconstructions<br>
   
  
- ***physicsframe.f90*** (module $\color{teal}\texttt{MOD\\_PhysicsFrame}$): Is where the initial and boundary conditions are defined, as well as the bathymetry function (and its derivative) and the two auxiliary subroutines for converting variables from conservative to primitive (subroutine ConsToPrim) and from primitive to conservative (subroutine PrimToCons). 
- ***mesh.f90*** (module $\color{teal}\texttt{MOD\\_Mesh}$): All the information related to the mesh (for 1D or 2D horizontal domains) is created in this file, namely, the nodes, barycenters, normal and tangent vectors at the cell edges, define the Gauss quadrature rule, and the quadrature points on each cell used to compute high order integrals in the finite volume scheme. 
- ***timediscretization.f90*** (module $\color{teal}\texttt{MOD\\_TimeDiscretization}$): Time integrators with one or several stages are defined in this module. On each method, each stage calls the subroutine FVSpaceIteration which performs the spatial approximation. 
- ***finitevolume.f90*** (module $\color{teal}\texttt{MOD\\_FiniteVolume}$): The pricipal subroutine in this module is FVSpaceIteration, which contains the spatial loops and contains the numerical methods used to approximate the conservative flux, non-conservative fluctuations and explicit approach of the source term. This module also calls the corresponding matrices of the system and source vector from the module MOD_MomentModels.
- ***reconstruction.f90*** (module $\color{teal}\texttt{MOD\\_Reconstruction}$): This module is dedicated to compute the high-order polynomial reconstructions at each cell edge of the mesh. A subroutine for slope limiters is also part of this module.
- ***shocksindicator.f90*** (module $\color{teal}\texttt{MOD\\_ShocksIndicator}$): Is used to identify the cells in which a shock occurs.
- ***momentmodels.f90*** (module $\color{teal}\texttt{MOD\\_MomentModels}$): The system matrix $`\boldsymbol{A}`$ and remaining matrix $`\boldsymbol{B}`$, as well as the conservative flux $`\boldsymbol{F}`$ and source vector $`\boldsymbol{S}`$ are defined in this module. 
- ***nonlinearsolver.f90*** (module $\color{teal}\texttt{MOD\\_NonlinearSolver}$): The Newton-Raphson solver is implemented in this module to solve the source term by means of an implicit scheme in an splitting procedure.
- ***output.f90*** (module $\color{teal}\texttt{MOD\\_Output}$): Dedicated to write the solution, mesh and bathymetry to ".dat" files.
- ***mytests.f90*** (module $\color{teal}\texttt{MOD\\_MyTests}$): This module is created for testing new subroutines on any of the modules of the program. It can be called from any of the subroutines or the main program. Module for developing purposes only.


## Running the program

After compilation, the executable program **main** is created in the folder *bin*, and to run a specific example, a folder related to the example needs to be created in *bin*, for instance, we create the folder *MyExample*. Then, the example folder needs to contain a file called **input.txt** which has a specific format for setting the parameters and flags of our example, described next. The example can be run at the bin folder as follows:

```console
user@pc:~/shallow-water-moments/bin$ ./main MyExample
```

The input file contains relevant information to run the specific example which can be modified without the need of compiling the program. At the input file, the user can set the model, numerical scheme and parameters, among others, see [**input.txt**](https://github.com/juliocareaga/shallow-water-moments/blob/main/bin/MyExample/input.txt). Please note that lines in the input file cannot be exchanged and real numbers need to have a decimal separator, while integers must not contain it. On each line the chosen number (single number, real or integer) needs to be typed once before the exclamation mark.

## Output data

The output data is created in the example folder. The solution is written in a sequence of ".dat" files with the following format **soln** followed by an integer which is the number of saved solution, and either **_1D** or **_2D** depending if the simulation is in one or two spatial dimensions, respectively. For instance, **soln3_1D.dat** corresponds to the third saved solution of a one dimensional model. The time points related to the sequence of saved solutions are saved in the file **timepoints.dat**. A file called **information.txt** containing relevant information of the computed simulation is also created at the example folder. Finally, when an exact solution is available, the $`L^1`$-error of each saved snapshot are created in a subfolder called *Error*.

## Authorship

This Fortran solver has been developed by **Julio Careaga** (https://github.com/juliocareaga/) with the collaboration of:
**Mirco Ciallella**, **Julian Koellermeier**, **Afroja Parvin** and **Rik Verbiest**. The main structure of the code, including the makefile and a number of subroutines in this program were partially based on the shallow water solver fv-solver-sw (https://github.com/jbnunezd/fv-solver-sw). 

## Papers:

To acknowledge the use of this software, cite the (current) paper:

$\color{blue}\texttt{(Current version)}$
J. Careaga, Q. Huang and J. Koellermeier. **A moment model of shallow granular flows with variable friction laws**, 
*arXiv preprint arXiv:2512.15332*, 2025, https://doi.org/10.48550/arXiv.2512.15332


