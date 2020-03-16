###########################################################
# Created: Mon Dec 17 15:38:05 2012
#
# Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
#
###########################################################


algs=("2-PR" "2-PSR" "2-PT" "2-PST" "2-PRT" "2-PSRT" "2-PRx" "2-PSRx" "2-PTx" "2-PSTx" "2-PRTx" "2-PSRTx")

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
        ninst=$(fat $j)
        ./prog -a $a -n $j -x 0 -s 0 -q $ninst -i $fname -v 1 >> $a.out
        mv $a.$j results-$a
    done
done
