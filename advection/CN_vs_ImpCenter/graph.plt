if(exist("n")==0 || n<0) n = n0
title(n) = sprintf("t = %.1f",n*5.0)
unset label
set label title(n) font 'Times,20' at 1.0, 1.0

plot "adv_imp_center.dat" index n ti "Center" w l, \
     "adv_CN.dat" index n ti "CrankNicolson" w l

n = n + dn
if (n < n1) reread
undefine n