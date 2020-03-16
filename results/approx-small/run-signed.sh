###########################################################
# Created: Mon Dec 17 15:38:05 2012
#
# Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
#
###########################################################


algs=("2-SPR" "2-SPSR" "2-SPRT" "2-SPSRT" "2-SPRx" "2-SPSRx" "2-SPRTx" "2-SPSRTx")

maxsize=${1:?"max_size path_instances"}
folder=${2:?"max_size path_instances"}

function fat() {
    local i=$1
    local f
    declare -i i
    declare -i f
    [ $i -le 2 ] && echo $i || {
        f=$(( i - 1));
        f=$(fat $f);
        f=$(( f * i ));
        echo $f;
    }
}

function expo() {
    set $1 $2 1
    while [ $2 -gt 0 ]; do
      set $1 $(($2-1)) $(($1*$3))
    done
    echo $3
}

if [ ! -f "prog" ]; then
    cd ../../
    make
    make clean
    mv prog results/approx-small
    cd results/approx-small
fi
    
for a in ${algs[@]}; do
    if [ ! -d "results-$a" ]; then
        mkdir results-$a
    fi
done

for ((j=2; j<=$maxsize; j++)); do
    for a in ${algs[@]}; do
        fname="$folder/$j.txt"
        ninst=$(( $(fat $j) * $(expo 2 $j) ))
        ./prog -a $a -n $j -x 0 -s 1 -q $ninst -i $fname -v 1 >> $a.out
        mv $a.$j results-$a
    done
done
