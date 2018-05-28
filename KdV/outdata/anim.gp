reset
set nokey
set xrange [0:30]
set yrange [-2:3]

set term gif animate
set output "u_anim.gif" 

n0 = 0
n1 = 200
dn = 1

load "graph.plt"