set terminal epslatex standalone color colortext 12
set output 'p01.tex'

set ytics nomirror
#set yrange [0:1400]
#set y2range [-30:40]
#set y2tics 10
#set ytics 10

s0 = 932.6173203
set lmargin 9
set rmargin 8
set x2tics 200
unset xtics
set multiplot layout 2,1 rowsfirst

# Set to fit with the plot
set size 1,0.2
set origin 0,0.8

#unset xrange 
set yrange [-1.3:1.3]
#set xrange [0:100]
unset ytics
set mxtics 5
set xzeroaxis 
set style fill solid
#unset ylabel
#unset y2label

unset key
#set lmargin 10
#set rmargin 8
set zeroaxis

# Elements on the to plot #
# 8 -> angle
# 6 -> k1l
# 7 -> k2l
# select, flag=twiss, column=NAME, KEYWORD, S, L, BETX, BETY, ALFX, ALFY,
# |MUX, MUY, DX, DPX, DY, DPY, ANGLE,K1L, K2L, K3L, K4L, envx, envy,k0l;
s = 3
l = 4
betx = 5
bety = 6
alfx = 7
alfy = 8
mux= 9
muy = 10
dx = 11
dpx = 12
dy = 13
dpy = 14
angle = 15
k1l = 16
k2l = 17
k3l = 18
k4l = 19
envx = 20
envy = 21
k0l = 22
set x2range [-1000:600]
set xrange [-1000:600]
unset x2tics

plot '< grep -i BEND ../mditrackall.tfs' u ($3-s0):($15 != 0 ? -0.5:1/0):($4) w boxes axis x2y1 lt 3 lw 1,\
     '< grep -i BEND ../mditrackall.tfs' u ($3-s0):($15 != 0 ?  0.5:1/0):($4) w boxes axis x2y1 lt 3 lw 1,\
     '< grep -i QUADRUPOLE ../mditrackall.tfs' u ($3-s0):($16/abs($16)):($4) w boxes axis x2y1 lt 1 lw 1,\
     '< grep -i SEXTUPOLE ../mditrackall.tfs' u ($3-s0):($7/abs($7)):($4+4) w boxes axis x1y1 lt 0 lw 1
#     '< grep -i KICKER ../mditrackall.tfs' u ($3-$4):(1):(0.2) w boxes axis x1y1 lt 4 lw 1



# plot '../tle_ae.tls' u 2:8 w lp ti 'yenvelope fw', '../tle_ae_rev.tls' u 2:8 w lp ti 'yenvelope bw'
unset yrange
set ytics nomirror
set y2tics
unset x2tics
set xtics 200
set size 1,0.9
set key

# Plot envelope
#plot '../tle_ae.tls' u 3:20 w lp ti 'henvelope fw', '../tle_ae_rev.tls' u (s0-$3):20 w lp ti 'henvelope bw'

# Plot \betax, \etax
set xlabel '$s$ [m]'
set ylabel '$\beta$ [m]'
set y2label '$\eta$ [m]'
set mytics 5
set my2tics 5
set y2range [-1:]

set key top left
#plot '../tle_ae.tls' u 3:5 w lp ti '$\beta_x$ fw', '../tle_ae_rev.tls' u (s0-$3):5 w lp ti '$\beta_x$ bw', '../tle_ae.tls' u 3:11 w lp ti '$\eta_x$ fw' axes x1y2, '../tle_ae_rev.tls' u (s0-$3):13 w lp ti '$\eta_x$ bw' axes x1y2

# # Plot \betay, \etay
# set xlabel '$s$ [m]'
# set ylabel '$\beta$ [m]'
# set y2label '$\eta$ [m]'
# set mytics 5
# set my2tics 5
# plot '../tle_ae.tls' u s:bety w lp ti '$\beta_y$ fw', '../tle_ae_rev.tls' u (s0-column(s)):bety w lp ti '$\beta_y$ bw', '../tle_ae.tls' u s:dy w lp ti '$\eta_y$ fw' axes x1y2, '../tle_ae_rev.tls' u (s0-column(s)):dy w lp ti '$\eta_y$ bw' axes x1y2

plot '../mditrackall.tfs' u ($3-s0):($5*1.0e-0) w lp ti '$\beta_x$', '' u ($3-s0):($6*1e0) w lp ti '$\beta_y$', '' u ($3-s0):11 w lp ti '$\eta_x$' axis x1y2
