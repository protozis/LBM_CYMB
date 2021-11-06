# Lattice Boltzmann Method - Cylindrical Moving Boundary

## Script dependences
These dependences are only used by Bash script, has no effect on the main simulation program.
- `time`: Linux built-in one, GNU version also works.
- `ffmpeg`: used by `video_maker`, pack .png images into .mp4 video.

## Build from source
Most of the C program is written in `C99` standard, therefore it should not need any extra lib. However, since some of the environment setups would be nasty for `clang` when you are compiling kernel program, I will recommend you to use `glibc` instead. As for the The OpenCL driver, it really dependent on the platform you have. You should check your OS instruction manual for which Opencl driver package you should install. In Archlinux they are:

**Runtime**
- Intel GPU: `intel-compute-runtime`
- Intel CPU: `intel-opencl-runtime<sup>AUR</sup>`
- Nvidia GPU: `opencl-nvidia`
- AMD GPU: `opencl-mesa`
- AND CPU: Not supported anymore.(?)

**Development**
- ICD loader: `ocl-icd`
- Headers: `opencl-headers`

**Tools**
- `clinfo`: good for finding all possible properties of the system.

```shell
$> cd src
$> make
$> make install
```

Be advice that the Opencl target version need to be defined in you host program as 
```C
#define CL_TARGET_OPENCL_VERSION 300
```
300 stands for Ver. 3.0.0

The atomic function support for 64-bits integer is required by the force calculation, need to check the device extension of `cl_khr_int64_atomics` of desire platform.

## Usage
Following steps can be altered by your own needs, feel free to play around with it.

### Choose a Resolution
To decide a good resolution for the simulation, you should consider with:
1. **The fluid phenomenon you are studying at.** The higher the resolution is, the more accurate simulation you get. However there dose not have a good calculation to tell you what resolution is enough for you, the only approach is performing a brunch of test run.
2. **The device you are using.** All devices will have their own `preferred work group size multiple` affected by the number of compute units and the size of cache. Check the value with `clinfo`. Most of the case, Intel CPU will go for `128` multiple, while GPU will go for `32` multiple. A simple approach is following the screen resolution, since it is how GPU is designed for. However it may not always be the best one for sure.

### Make an experiment setup
Inside an experiment setup directory we have:
```
example/
├── a.bc
├── a.nd
├── a.pd
├── data
│   ├── 0.mp4
│   ├── 1.mp4
│   └── 2.mp4
├── debug
├── default.conf
├── log
└── plot.p
```

|Name|Description|
|-|-|
|a.bc|Boundary condition file|
|a.nd|Number Density file|
|a.pd|Platform/Device file|
|data|Result videos or ND files|
|debug|Used to collect desire parameters|
|log|Running conditions and status|
|plot.p|gnuplot script|

#### ND (Number Density) file
Number density matrix represent particle densities in different velocities and positions. For an D2Q9 ND file with 480x270 in size, is formatted as
```
480 270 9
[m0q0] [m0q1] [m0q2] ...
[m1q0] [m1q1] [m1q2] ...
[m2q0] [m2q1] [m2q2] ...
...
```

Use `nd_gen` to make an initial ND file with following options:
```
Usage: nd_gen [OPTION]
Otions:
	-x #, Numer of column
	-y #, Number of row
	-d #, Density (Default: 0.3)
	-i #, Ux (Default: 0.12)
	-j #, Uy (Default: 0)
	-h , Print this help page
Example:
	nd_gen -x 100 -y 30 -d 0.3 -i 0.1 -j 0 > exp_sets/[EXP_NAME]

```  

#### BC (Boundary Condition) file
Boundary condition consist two data types:
1. BCV: Unified fluid that surround the simulation area.
2. CY: Cylinder object.

For example:
```
bc_no 2
bc_nq 2
BCV {
	dnt 1
	ux 0.0
	uy 0.0
}
CY {
	spring 0.1
	damp 0
	mass 1000
	rad 20
	force 0 0
	acc 0 0
	vel 0 0
	dsp 10 0
	pos 250 70
}
CY {
	spring 0.1
	damp 0
	mass 1000
	rad 20
	force 0 0
	acc 0 0
	vel 0 0
	dsp 10 0
	pos 250 150
}
```
|Parameter|Description|
|-|-|
|bc_no|Number of objects|
|bc_nq|Kinematic dimension of objects, usually 2 in 2D simulation.|
|dnt|Density|
|ux|Macro velocity in x axis|
|uy|Macro velocity in y axis|
|spring|Spring constant of the cylinder, in lattice unit|
|force|Initial value of applied force of the cylinder, in lattice unit|
|damp|Damping constant of the cylinder, in lattice unit|
|mass|Mass of the cylinder, in lattice unit|
|rad|Radius of the cylinder, in lattice unit|
|force|Initial value of applied force of the cylinder, in lattice unit|
|acc|Initial value of applied force of the cylinder, in lattice unit|
|vel|Initial value of velocities of the cylinder, in lattice unit|
|dsp|Initial value of displacements of the cylinder, in lattice unit|
|pos|Initial value of positions of the cylinder, in lattice unit|

#### PD (Platform/Device) file

> **Work-group & Work-item** 
> A work-group is processed by a single compute unit in the device. For a CPU device it prefer a larger work-group size, while a GPU works opposite. The amount of work-item in a work-group has its limitation, check *Max work group size* for the value. For multi-dimension work-group, the total amount of work-item cannot exceed *Max work group size* too, and the work-item in each dimension will have their own limitation. Check *Max work item sizes* for the value. However, I will recommend to choose the work-group size with automatic tools `speed_test`.

Use the script `speed_test` to choose the best environment setup. If the chosen work-group exceed the max allowed work-group size for the device, the result will not be printed. Example: 

```shell
$> ./speed_test -i exp_sets/example
#Syntax: pdFileName[platform device localSize1 localSize2] realTime(secs)
#
Working on...completed!
p2d0_120.pd[2/0/16/9]                             	7.90
p0d0_120.pd[0/0/16/9]                             	8.11
p1d0_24.pd[1/0/80/45]                             	8.28
p1d0_40.pd[1/0/48/27]                             	8.43
p1d0_60.pd[1/0/32/18]                             	8.52
p1d0_120.pd[1/0/16/9]                             	8.55
p1d0_20.pd[1/0/96/54]                             	8.60
p1d0_30.pd[1/0/64/36]                             	8.70
p2d0_60.pd[2/0/32/18]                             	8.75
```
You can see that Nvidia GPU with 32x18 work-group size is the fastest one. Write the result to your experiment setup with:
```
$> cp pd_auto/p2d0_120.pd exp_sets/[EXP_NAME]/a.pd
``` 


### Reference
#### Scripts
|Name|Purpose|
|---|---|
|test_run|LBM simulation wrapper (with '-' argument)|
|scheme_run|LBM simulation wrapper (with configuration file)|
|speed_test|PD file speed tester|
|video_maker|Generate mp4 from png|
|maker_bc|Generate BC file from PNG_BC|
|plot_obj|Plot object force diagram|
|run_plot_sp|Plot speed profile for each column|
|run_plot_sum|Plot flow rate for each cross section|
|viewer_nd|Image view a ND file|
|viewer_bc|Image view a BC file|
|viewer_sp|Plot speed profile|
|viewer_sum|Plot flow rate|
|viewer+total|Plot total particles|
|pd_gen_auto|Auto generating pd files|
#### Binaries
|Name|Purpose|
|---|---|
|simulator|Opencl LBM simulation with '-' argument|
|launcher|Opencl LBM simulation with configuration file|
|nd_gen|Generate initial ND file|
|nd_get_par|Get macro-scale properties|
|nd_get_sum|Get total particle|
|bc_gen_box|Generate BC file with a box boundary|
|bc_gen_tunnel|Generate BC file with a wind tunnel boundary|
|bc_gen_field|Generate BC file with a uniform flow|
|bc2ppm|Convert BC file to ppm image|
|ppm2bc|Convert PPM image to BC file|
#### Files
|Name|.|Purpose|Format|
|:---:|:---:|---|----|
|ND|.nd|Number density matrix|[nx] [ny] [nq]\n[m]; addr = q+(x+y*nx)*nq; (double #)number_density|
|BC|.bc|Boundary matrix|[nx+2] [ny+2]\n[m]; addr = x+y*nx; (char [0-9])Object/ (char '>')input/ (char '#')output/ (char '+')fluid ;|
|PD|.pd|Platform/Device profile|[Platform] [Device] [Workgroup dimension 1] [workgroup dimension 2]|
|SVG_BC|.svg|Boundary image|(#000000)fluid/(#ffff09)border/ (#ffff01)obj_1/ (#ffff02)obj_2...|
|OBJ||Object force|[loop] [Fx] [Fy]|
