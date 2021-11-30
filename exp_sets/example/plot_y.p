#! /usr/bin/gnuplot

set grid
plot 'debug' u 2:6 w lp, '' u 2:8 w lp, '' u 2:10 w lp
pause mouse close
