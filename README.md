# Lattice Boltzmann Method - Cylindrical Moving Boundary

![3cir](img/3cir.gif)

Before you start using this simulator, it is highly recommended to read the instruction I made: **[What's the physics of this LBM simulation?](physics.md)**.

## Table of contents
1. [Dependences and Build process](#Dependences-and-Build-process)
2. [Start a Simulation Step by Step](#Start-a-Simulation-Step-by-Step)
3. [Design and Analyze an Experiment](#Design-and-Analyze-an-Experiment)

## Dependences and Build process
### C Binaries
Most C programs are written in `C99` standard, therefore no extra libs needed. However, since some of the environment setups would be nasty for `clang` when you are compiling OpenCL kernel program, I will recommend you to use `glibc` instead. As for the The OpenCL driver, it really dependent on the platform you have. You should check your OS instruction manual for the driver packages you need to install. In Archlinux they are:

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
- `clinfo`: good for monitoring all possible platform properties of the system.

You can build and install the binaries with: 
```shell
$> cd LBM_CYMB/src
$> make
$> make install
```
`make install` will copy all binaries and `.cl` kernel source file into `LBM_CYMB/bin`. Clean all binaries with `make clean` if you want a fresh make.

Be advice that the Opencl target version need to be defined in you host program as 
```C
#define CL_TARGET_OPENCL_VERSION 300
```
while `300` stands for the version `3.0.0`. Also the atomic function support for 64-bits integer is required by the force calculation in `simulate_ocl.cl`, you will need to check if the device extension of `cl_khr_int64_atomics` is available for your desire platform.

> [Note] If you want to make MP4 video in pipeline, `ffmpeg` must be installed. 

### Scripts
These dependences are only used by shell script.
- `time`: Linux built-in one, GNU version also works.
- `gnuplot`: Used by `plot.p` for visualize data.

## Start a Simulation Step by Step
Following steps can be altered by your own needs, feel free to play around with it.

### Choose a Resolution
To decide a good resolution for the simulation, you should consider with:
1. **The fluid phenomenon you are studying at.** 

	The higher the resolution is, the more accurate simulation you get. However there dose not have a good calculation to tell you what resolution is enough for you, the only approach is benchmarking.

2. **The device you are using.**

	All devices will have their own `preferred work group size multiple` affected by the number of compute units and the size of cache. Check the value with `clinfo`. Most of the case, Intel CPU will go for `128` multiple, while GPU will go for `32` multiple. A simple approach is following the screen resolution, since it is how GPU is designed for. However it may not always be the best one for sure.

### Make an Experiment Setup
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

To perform an experiment, you will need to tweak these files:
|Name|Description|
|-|-|
|a.bc|Boundary condition file|
|a.nd|Number Density file|
|a.pd|Platform/Device file|
|default.conf|Main configuration file|

And following files are used to collect and analyze data.
|Name|Description|
|-|-|
|debug|Used to collect desire parameters|
|data|Result videos or ND files|
|log|Running conditions and status|
|plot.p|gnuplot script|


#### ND (Number Density) file
Number density matrix describe particle densities in different velocities and positions. For an D2Q9 ND file with 480x270 in size, is formatted as
```
480 270 9
[m0q0] [m0q1] [m0q2] ...
[m1q0] [m1q1] [m1q2] ...
[m2q0] [m2q1] [m2q2] ...
...
```
![ND](img/nd.png)

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
	nd_gen -x 100 -y 30 -d 0.3 -i 0.1 -j 0 > exp_sets/example/a.nd

```  

#### BC (Boundary Condition) file
Boundary condition include following classes:
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

- **Platform & Device**
	
	A platform is an specific OpenCL implementation, e.g. Intel, AMD or Nvida CUDA. A device is the actual processor that perform the calculation, like `Intel(R) Core(TM) i7-8565U CPU @ 1.80GHz`, `NVIDIA GeForce MX150`.

- **Work-group & Work-item** 

	A work-group is processed by a single compute unit in the device. For a CPU device it prefer a larger work-group size, while a GPU works opposite. The amount of work-item in a work-group has its limitation, check *Max work group size* for the value. For multi-dimension work-group, the total amount of work-item cannot exceed *Max work group size* too, and the work-item in each dimension will have their own limitation. Check *Max work item sizes* for the value.

Use the script `speed_test` to choose the best environment setup. If the chosen work-group exceed the max allowed work-group size for the device, the result will not be printed. Example: 

```shell
$> clinfo -l
Platform #0: Intel(R) OpenCL HD Graphics
 `-- Device #0: Intel(R) UHD Graphics 620 [0x3ea0]
Platform #1: Intel(R) OpenCL
 `-- Device #0: Intel(R) Core(TM) i7-8565U CPU @ 1.80GHz
Platform #2: NVIDIA CUDA
 `-- Device #0: NVIDIA GeForce MX150
$> ./speed_test exp_sets/example
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
You can see that Nvidia GPU with 16x9 work-group size is the fastest one. Write the result to your experiment setup with:
```
$> cp pd_auto/p2d0_120.pd exp_sets/example/a.pd
``` 

#### .conf file
Following parameters are included in a configuration file:
|Parameter|Description|Default|
|-|-|-|
|LOOP|Iteration of `SKP`.|1|
|SKP|Iteration of fluid propagate by time step.|1|
|CF|Collision frequency in Lattice unit|1|
|CS|Speed of sound in SI unit.|340|
|CL|Dimensional value of length|1|
|CD|Dimensional value of Density|1|
|MA|Mach number(Dimensionless)|0.1|
|IS_MP4|Make mp4 videos or not.|0|
|IS_SAVE_DATA|Save final ND matrix or not.|0|
|IS_FILE_OUTPUT|Save ND matrix for every loops or not.|0|
|ND_FILE|ND filename|NULL|
|BC_FILE|BC filename|NULL|
|PD_FILE|PD filename|NULL|
|OUTPUT_DIR|Output filename|NULL|
|PROGRAM_FILE|OpenCL source filename|NULL|
|REFUEL_RTO|Refuel ratio|0.8|
|EAT_RTO|Eat ratio|0.01|
|LOG_FILE|Log filename|NULL|
|DEBUG_FILE|Debug filename|NULL|
|IS_LOG_PRINT|Print log to stdout or not.|1|
|IS_PROGRESS_PRINT|Print progress or not|1|
|PL_MAX_D|Maximum density value for the jetcolormap ploting|0.5|
|PL_MAX_UX|Maximum x velocity value for the jetcolormap ploting|0.1|
|PL_MAX_DUY|Maximum y velocity value for the jetcolormap ploting|0.1|

### Perform an experiment
Run `simulator` with desire experiment setup:
```shell
$> ./simulator exp_sets/1920x1080
``` 
If multiple experiments are going to be launched at once, use `schedule` plus a plain text list file of experiment setups directory name, like:
```shell
$> cat exp_grp1
exp_sets/480x270_skp10
exp_sets/480x270_skp20
exp_sets/480x270_skp30
exp_sets/480x270_skp40
exp_sets/480x270_skp50
```
And launch with:
```shell
$> ./schedule exp_grp1
```

## Design and Analyze an Experiment
Here I will show you a working example for designing and analyzing an experiments. Today I am curious about the behaviors of two cylinders listed alongside the flow direction, 
with different offsets between them as shown below:

![exp_2c-1](img/exp_2c.gif)

![exp_2c-2](img/exp_2c.png)

