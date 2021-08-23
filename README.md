## Lattice Boltzmann Method Simulator
![icon](icon.gif)
### Patch notes
#### 2020-06-09
1. Removed pthread LBM simulation.
>The performance is far too weak compare with OpenCL.
2. Schedule simulation with configuration file.
>Use `.conf` for rapid launch.
3. Streaming produced MP4 for image view.
>Reducing storage usage to limit.
4. Atomic function(C11) included.
>Conquered data racing problem in opencl parallelizing.
5. Moving boundary support.
>That's it, finally.

### Branches
#### Master
Single device with 2 dimensions work-groups. Since the processing time for 1920x1080 number density matrix for 200 loop times is less then 12 secs in all devices on my laptop, this branch perform different fluid parameter setups on each device independently.

#### 1-D_with_eq_scaling
This kernel working in 1 dimension work-group. In some configurations in CPU device it will be slightly faster then 2D work-group, but will be slower for all configurations in GPU device.

#### 1PlatformNdevice(Working in progress)
The context in this branch content more then 1 device in a single platform, e.g. server motherboard with multiple CPU on it. Since there is only 1 RAM for every devices, higher processing speed is expected. However most PC will only have 1 CPU, and all GPU card are treated as 1 device with 1 internal RAM, this branch will has less primary then master.

#### NplatformNdevice(Not initiated yet)
Divide the number density matrix for each device on different platform. However due to the calculations of macro-scale parameters and discretized distribution function, the edge of result number density need to be exchanged between each devices on every time steps. The data will need to be transferred between motherboard RAM and GPU devices RAMs constantly, spending plenty of times and computing resources.

### For every scripts and binaries...
All executable files have manual page. Print it with `-h` option or simply left no argument (not everyone do this trick). For example:
```shell
$> ./test_run -h
Usage: test_run [OPTION]

Option:
        -i [ND_FILE], input ND file
        -b [BC_FILE], input BC file
        -o [DIR], output directory
        -n #, loop number
        -p #, plot png data
        -v #, plot data to mp4 video
        -t #, skip ts (default 100)
        -c #, collision (default 1)
        -l #, lattice speed (default 1)
        -k #, pd file path
        -h , Print this help page

```
### Script dependences
This dependences are only used by Bash script, has no effect on the main simulation program.
- `time`: Linux built-in one, GNU version also works.
- `feh`: Image viewer.
- `ffmpeg`: used by `video_maker`, pack .png images into .mp4 video.
- `inkscape`: a good tool for generating svg files.

### Build from source
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

For pthread version, it should also come with OS since it is included in standrad POSIX. However this version is only used for debug purpose, the speed is much much slower then Opencl version, so **avoid using it...**

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

### Usage
Following steps can be altered by your own needs, feel free to play around with it.

#### Choose a Resolution
To decide a good resolution for the simulation, you should consider with:
1. **The fluid phenomenon you are studying at.** The higher the resolution is, the more accurate simulation you get. However there dose not have a good calculation to tell you what resolution is enough for you, the only approach is performing a brunch of test run.
2. **The device you are using.** All devices will have their own `preferred work group size multiple` affected by the number of compute units and the size of cache. Check the value with `clinfo`. Most of the case, Intel CPU will go for `128` multiple, while GPU will go for `32` multiple. A simple approach is following the screen resolution, since it is how GPU is designed for. However it may not always be the best one for sure.

#### Generate ND and BC files
1. Use `nd_gen` to make an initial ND file.
2. Use `inkscape` to plot a boundary in svg format with the same resolution. Check [References](#files) for detail.
3. Make the BC file from svg file with `maker_bc`

#### Choose a work-group size
A work-group is processed by a single compute unit in the device. For a CPU device it prefer a larger work-group size, while a GPU works opposite. The amount of work-item in a work-group has its limitation, check *Max work group size* for the value. For multi-dimension work-group, the total amount of work-item cannot exceed *Max work group size* too, and the work-item in each dimension will have their own limitation. Check *Max work item sizes* for the value. However, I will recommend to choose the work-group size with automatic tools `speed_test` introduced in next section.

#### Generate PD file & Speed Test
1. Use `pd_file/pd_gen` to generate all possible PD file (Platform/Device file). The global size most be divide perfectly by local size, while the generator will handle it automatically. For example if you are working on 1920x1080 resolution, generate PD files with:
```shell
$> ./pd_gen -x 1920 -y 1080 -p 0 -d 0 -n device_p0d0 -m 1080
```
Or, you can auto generate PD file with `pd_gen_auto`.

2. Use the script `speed_test` to choose the best environment setup. If the chosen work-group exceed the max allowed work-group size for the device, the result will not be printed. Example: 

Auto generate PD file with `pd_auto_gen`
```shell
$> ./speed_test -i nd_files/1920x1080.nd -b bc_files/1920x1080_2c.bc -o data -t 20
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
Pre-generated PD files
```shell
$> ./speed_test -i nd_files/1920x1080.nd -b bc_files/1920x1080.bc -o data -t 50 -k pd_file/
#Syntax: pdFileName[platform device localSize1 localSize2] realTime(secs)
#
Working on...completed!
nv_gpu_60.pd[2/0/32/18]                           	9.46
nv_gpu_120.pd[2/0/16/9]                           	9.56
intel_gpu_120.pd[0/0/16/9]                        	9.76
intel_cpu_24.pd[1/0/80/45]                        	10.07
intel_cpu_40.pd[1/0/48/27]                        	10.11
intel_cpu_120.pd[1/0/16/9]                        	10.16
intel_cpu_30.pd[1/0/64/36]                        	10.20
intel_cpu_60.pd[1/0/32/18]                        	10.26
```

You can see that Nvidia GPU with 32x18 work-group size is the fastest one. Next time you can perform the test run with the PD file `nv_gpu_60.pd` directly.

#### Test Run
If the number density matrix and the object forces are the only data you need for the simulation result (Usually they are), initiate the test run with the best PD file you got from `speed_test` with:
```shell
$> ./test_run -i nd_file/1920x1080.nd -b bc_file/1920x1080.bc -o data -n 100 -t 10 -k pd_file/nv_gpu_60.pd
```
How to choose a good loop-time/skip-time combination? A higher loop-time will give you smoother simulation, as for a higher skip-time will be faster in process because of the reducing of out stream.

#### Schemdule Run
After all environment setups and parameters are decided, you can used `scheme_run` with a configuration file to perform your simulation in a more rapid way. For example, `example.conf` can be found in `configs` :
```
LOOP	10
SKP	10
CF	1
SL	1
IS_MP4	1
IS_SAVE_DATA	0
IS_FILE_OUTPUT	0
OBJ_MAX	10
PLOT_MAX	0.5
FC_OFS	1000000
FC_RIST	0
FC_DAMP	0
FC_MASS	0.05
FC_DT	0.001
FC_NQ	2
ND_FILE	nd_files/1920x1080.nd
BC_FILE	bc_files/1920x1080_2c.bc
PD_FILE	pd_file/1920x1080/nv_gpu_120.pd
DIR_NAME data
```
And now you can launch your simulation with:
```shell
$> ./scheme_run example.conf
```

#### Analyzing Data
1. **(NOT RECOMMEND!)** For image viewing of the jet-color field, plot the result with `-p` option. You can also combine the images to videos with `video_maker` like this:
```shell
$> ./video_maker data/plot/combine output.mp4
```
2. Since plotting all image in high resolution will cost huge storage space (>100G), it is recommended to stream-produce MP4 video without saving raw image. You can adding `-v` when using `simulator` in `test_run`, or adding `IS_MP4 1` when using `launcher` in `scheme_run`.

3. The result of object forces can be shown with `plot_obj`. If you have 4 objects you can plot the force in x axes with:
```shell
$> ./plot_obj -i data -d 2 -n 4
```

4. Two analyzing tools `run_plot_sum` and `run_plot_sp` are used in *poiseuille flow* test, can plot speed profile of each column and fluid rate in a set of cross sections. `viewer_total` is a simple script that return the total particle numbers inside the simulation space in each time step.

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
