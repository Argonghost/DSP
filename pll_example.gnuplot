reset                                   # reset
set size ratio 0.2                      # set relative size of plots
set xlabel 'Sample Index'               # set x-axis label for all plots
set grid xtics ytics                    # grid: enable both x and y lines
set grid lt 1 lc rgb '#cccccc' lw 1     # grid: thin gray lines
set multiplot layout 3,1 scale 1.0,1.0  # set three plots for this figure

# real
set ylabel 'real'                       # set y-axis label
set yrange [-1.2:1.2]                   # set plot range
plot 'pll_example.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#999999' notitle ,\
     'pll_example.dat' using 1:4 with lines lt 1 lw 2 lc rgb '#004080' notitle

# imag
set ylabel 'imag'                       # set y-axis label
set yrange [-1.2:1.2]                   # set plot range
plot 'pll_example.dat' using 1:3 with lines lt 1 lw 2 lc rgb '#999999' notitle ,\
     'pll_example.dat' using 1:5 with lines lt 1 lw 2 lc rgb '#008040' notitle

# phase error
set ylabel 'phase error'                # set y-axis label
set yrange [-3.2:3.2]                   # set plot range
plot 'pll_example.dat' using 1:6 with lines lt 1 lw 3 lc rgb '#800000' notitle

unset multiplot