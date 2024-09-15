reset                                   # reset
set size ratio 0.2                      # set relative size of plots
set xlabel 'Sample Index'               # set x-axis label for all plots
set grid xtics ytics                    # grid: enable both x and y lines
set grid lt 1 lc rgb '#cccccc' lw 1     # grid: thin gray lines
set multiplot layout 3,1 scale 1.0,3.0 

# real
set ylabel 'Nav bits'                       # set y-axis label
set yrange [-2:2]                   # set plot range
plot 'gps_data_output.txt' using 1:3 with lines lt 1 lw 3 lc rgb '#999999' notitle

# real
set ylabel 'CA Code'                       # set y-axis label
set yrange [-2:2]   
set xrange [0:300]                # set plot range
plot 'gps_data_output.txt' using 1:4 with lines lt 1 lw 2 lc rgb '#004080' notitle 

# complex
set ylabel 'GPS L1 Signal'                       # set y-axis label
set yrange [-2:2] 
set xrange [0:300]                   # set plot range
plot 'gps_data_output.txt' using 1:5 with lines lt 1 lw 2 lc rgb '#800000' notitle 

unset multiplot