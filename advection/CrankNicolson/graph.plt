if(exist("n")==0 || n<0) n = n0
title(n) = sprintf("t = %.1f",n*1.0)
unset label
set label title(n) font 'Times,20' at 1.0, 1.3

plot "adv.dat" index n w l

n = n + dn
if (n < n1) reread
undefine n