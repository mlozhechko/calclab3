#!/bin/bash

function tab1 {
    ./calc4 -precision double -grid uniform -fid $1 -report tab1_lag >> $2
}

function tab2lagC {
    ./calc4 -precision double -grid uniform -fid $1 -report tab2_lag_c >> $2
}

function tab2lagU {
    ./calc4 -precision double -grid uniform -fid $1 -report tab2_lag_u >> $2
}

function tab2spline {
    ./calc4 -precision double -grid uniform -fid $1 -report tab2_spline >> $2
}

function drawLag {
    for i in {2..7}
    do
        # echo "$i"
        ./calc4 -precision double -grid uniform -fid $1 -n $((2**i)) -method lagrange > /dev/null
        ./calc4 -precision double -grid chebyshev -fid $1 -n $((2**i)) -method lagrange > /dev/null
    done
}


function drawSpline {
    for i in {2..7}
    do
        # echo "$i"
        ./calc4 -precision double -grid uniform -fid $1 -n $((2**i)) -method spline > /dev/null
    done
}

echo "2 question"
tab1 const report/const.txt
drawLag const
tab1 linear report/linear.txt
drawLag linear

echo "3 question"
drawLag runge
drawSpline runge

echo "4 question"
drawSpline linear
drawSpline quad
tab2spline linear report/linear_spline.txt
tab2spline quad report/quad_spline.txt

echo "5 question"
drawSpline sin1
drawSpline sin125
drawLag sin1

tab2spline sin125 report/sin125_tab2_spline.txt
tab2lagC sin1 report/sin1_tab2_cheb.txt
tab2lagU sin1 report/sin1_tab2_unif.txt
tab2spline sin1 report/sin1_tab2_spline.txt

echo "6 question"
drawSpline target
drawLag target
tab2lagC target report/target_tab2_cheb.txt
tab2lagU target report/target_tab2_unif.txt
tab2spline target report/target_tab2_spline.txt