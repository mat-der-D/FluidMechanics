reset
set nokey
set xrange [0:50]
set yrange [-1.2:1.2]

set term gif animate
set output "adv_diff.gif" 

n0 = 0
n1 = 1000
dn = 1

load "graph.plt"