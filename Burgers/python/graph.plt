if(exist("n")==0 || n<0) n = n0
title(n) = sprintf("t = %.2f",n*0.01)
unset label
set label title(n) font 'Times,20' at 0.5, 1.8

plot "Burgers.dat" index n w l

n = n + dn
if (n < n1) reread
undefine n