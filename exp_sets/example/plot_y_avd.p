#! /usr/bin/gnuplot

set grid
plot 'debug' u 2:6 w linespoints, '' u 2:6 w linespoints, '' u 2:6 w linespoints
pause mouse close
