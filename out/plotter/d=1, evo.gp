#etaAMY_Nf0(g) = 27.126/(g**4*log(2.765/g))
#etaAMY_Nf2(g) = 86.47/(g**4*log(2.954/(sqrt(4/3)*g)))
#etaAMY_Nf3(g) = 106.66/(g**4*log(2.957/(sqrt(3/2)*g)))
#etaAMY_Nf6(g) = 147.63/(g**4*log(2.940/(sqrt(2)*g)))
eDen(t,c2s) = (1./t)**(1.+c2s)

set xl 'tau'
set yl "z"

set zr [0:5]

set grid
set key l b

set tit "d=1 hydro, Bj initial cond."
splot "../data/Evo.dat" u 1:2:3 with lines
pause -1

