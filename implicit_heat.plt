set colorbox user horizontal origin 0.32,0.385 size 0.18,0.035 front
set cbrange [236:1296]
set cbtics ('236 Hz' 236,'1296 Hz' 1296) offset 0,0.5 scale 0
set palette defined (\
    1  '#0025ad', \
    2  '#0042ad', \
    3  '#0060ad', \
    4  '#007cad', \
    5  '#0099ad', \
    6  '#00ada4', \
    7  '#00ad88', \
    8  '#00ad6b', \
    9  '#00ad4e', \
    10 '#00ad31', \
    11 '#00ad14', \
    12 '#09ad00' \
    )
plot for [n=1:6] "data.dat" u 1:(column(n)*1000) w lines ls n