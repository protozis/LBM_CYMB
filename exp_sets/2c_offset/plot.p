#! /usr/bin/gnuplot

var = 9

set grid
set title 'TITLE'
plot '< grep ^0 debug' u 2:var w lp, '< grep ^1 debug' u 2:var w lp
pause mouse close
