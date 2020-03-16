###########################################################
# Created: Qua 26 Ago 2015 16:52:05 BRT
#
# Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
#
###########################################################
 
algs=("WPRm" "WPRg" "WPR" "WPSRg" "WPSR" "WPTg" "WPT" "WPSTg" "WPST" "WPRTg" "WPRT" "WPSRTg" "WPSRT")

maxsize=${1:?"max_size path_instances alpha"}
folder=${2:?"max_size path_instances alpha"}
alpha=${3:?"max_size path_instances alpha"}

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

if [ ! -f "prog" ]; then
    cd ../../
    make
    make clean
    mv prog results/approx-small
    cd results/approx-small
fi

for a in ${algs[@]}; do
    if [ ! -d "results-${a}_${alpha}" ]; then
        mkdir results-${a}_${alpha}
    fi
done

for ((j=2; j<=$maxsize; j++)); do
    for a in ${algs[@]}; do
        fname="$folder/$j.txt"
        ninst=$(fat $j)
        ./prog -a $a -n $j -x $alpha -s 0 -q $ninst -i $fname -v 1 >> ${a}_${alpha}.out
        mv ${a}*.$j results-${a}_${alpha}
    done
done
