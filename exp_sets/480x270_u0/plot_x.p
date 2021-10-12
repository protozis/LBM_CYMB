#! /usr/bin/gnuplot

set grid
plot 'debug' u 2:3 w linespoints,'' u 2:5 w linespoints, '' u 2:7 w linespoints, '' u 2:9 w linespoints
pause mouse close
