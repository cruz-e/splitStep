 set term pslatex color solid

 set output  'fig.tex'
 set format "\\LARGE %\g"


 set multiplot
 set origin 0,.9
 set size 1.2,1.5
 set border lw 2

 set xlabel  "{\\LARGE \\bf $z$}" #  0,-1
 set ylabel  "{\\LARGE \\bf $T$}" # 0,1
 set zlabel offset 5,7 "{\\LARGE \\bf Power}"


set xrange[-10:10]
 set yrange[-10:10]
 set zrange[-0.1:1]

 set mxtics 2

 set mytics 2

 unset key 


 splot 'standard-fiber_power' u 1:2:5 w l lt 1 lw 4

