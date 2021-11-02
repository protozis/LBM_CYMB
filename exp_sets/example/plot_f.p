#! /usr/bin/gnuplot

set grid
plot 'debug' u 2:3 w lp, '' u 2:4 w lp
pause mouse close
