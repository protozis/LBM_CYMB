#! /usr/bin/gnuplot

set grid
plot 'debug' u 2:3 w line,'' u 2:5 w line, '' u 2:7 w line, '' u 2:9 w line
pause mouse close
