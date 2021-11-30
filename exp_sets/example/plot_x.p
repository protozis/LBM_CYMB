#! /usr/bin/gnuplot

set grid
plot 'debug' u 2:5 w lp, '' u 2:7 w lp, '' u 2:9 w lp
pause mouse close
