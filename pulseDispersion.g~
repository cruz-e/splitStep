 set term pslatex color solid
# set term x11
# OUTPUT
 set output  'fig.tex'
 set format "\\LARGE %\g"

# 
 set multiplot
 set origin 0,.9
 set size 1.2,1.5
 set border lw 2
#  
# set title "{\\Large \\bf Dispersi\\'on}"
 set xlabel  "{\\LARGE \\bf $z$}" #  0,-1
 set ylabel  "{\\LARGE \\bf $T$}" # 0,1
 set zlabel offset 5,7 "{\\LARGE \\bf Power}"
# 
 set auto
# set xrange[0:3]
 set yrange[-10:10]
 set zrange[-0.1:1]
# set xtics 0,1,3
 set mxtics 2
# set ytics -4,2,4
 set mytics 2
# set ztics 0,.5,1
# set grid
 unset key  # 1.5,-1.3

# plot 'standard-fiber_power' u 2:5 w l lt 1 lw 4


 splot 'standard-fiber_power' u 1:2:5 w l lt 1 lw 4

