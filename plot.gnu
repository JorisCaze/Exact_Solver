set xlabel("Abscissa x [m]")

set title("Density [kg/m3]")
plot "output.txt" u 1:2 w lp
pause(-1)

set title("Speed [m/s]")
plot "output.txt" u 1:3 w lp
pause(-1)

set title("Pressure [Pa]")
plot "output.txt" u 1:4 w lp
pause(-1)

set title("Mach [-]")
plot "output.txt" u 1:5 w lp
pause(-1)
