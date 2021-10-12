#! /usr/bin/gnuplot

set grid
plot 'debug' u 2:4 w line,'' u 2:6 w line, '' u 2:8 w line, '' u 2:10 w line
pause mouse close
