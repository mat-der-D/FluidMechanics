reset
set nokey
set xrange [0:50]
set yrange [-0.5:1.5]
# set yrange [-50:3000]

set term gif animate
set output "adv.gif" 

n0 = 0
n1 = 100
dn = 1

load "graph.plt"