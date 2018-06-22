reset
set nokey
set xrange [0:10]
set yrange [-1.0:2.0]

set term gif animate
set output "Burgers.gif" 

n0 = 0
n1 = 100
dn = 1

load "graph.plt"