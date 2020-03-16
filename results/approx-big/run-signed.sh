###########################################################
# Created: Mon Dec 17 15:38:05 2012
#
# Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
#
###########################################################

algs=("2-SPR" "2-SPSR" "2-SPRT" "2-SPSRT")

minsize=${1:?"min_size max_size path_instances nof_perms_per_file"}
maxsize=${2:?"min_size max_size path_instances nof_perms_per_file"}
folder=${3:?"min_size max_size path_instances nof_perms_per_file"}
ninst=${4:?"min_size max_size path_instances nof_perms_per_file"}

if [ ! -f "prog" ]; then
    cd ../../
    make
    make clean
    mv prog results/approx-big
    cd results/approx-big
fi

for a in ${algs[@]}; do
    if [ ! -d "results-$a" ]; then
        mkdir results-$a
    fi
done

for ((j=$minsize; j<=$maxsize; j+=5)); do
    for a in ${algs[@]}; do
        fname="$folder/$j.txt"
        ./prog -a $a -n $j -x 0 -s 1 -q $ninst -i $fname -v 1 >> $a.out
        mv $a.$j results-$a
    done
done

