Min = 0.0 # where binning starts
Max = 30.0 # where binning ends
n = 30 # the number of bins
binwidth = (Max-Min)/n # binwidth; evaluates to 1.0

bin(x) = binwidth*(floor((x-Min)/binwidth)+0.5) + Min

set boxwidth binwidth

set style histogram rowstacked gap 0
set style fill solid 0.5 border lt -1

p 'data.txt' u (bin($1)):(1.0) smooth freq w boxes lc rgb"blue" notitle

