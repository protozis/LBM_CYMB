#! /usr/bin/gnuplot

set grid
plot 'debug' u 2:4 w linespoints,'' u 2:6 w linespoints, '' u 2:8 w linespoints, '' u 2:10 w linespoints
pause mouse close
