if(exist("n")==0 || n<0) n = n0
title(n) = sprintf("t = %.1f",n*0.1)
unset label
set label title(n) font 'Times,20' at 0.2, 2.7

plot "u_anim.dat" index n w l

n = n + dn
if (n < n1) reread
undefine n