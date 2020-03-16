###########################################################
# Created: Mon Dec 17 15:38:05 2012
#
# Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
#
###########################################################

minsize=${1:?"min_size max_size path_instances nof_perms_per_file algorithm signal"}
maxsize=${2:?"min_size max_size path_instances nof_perms_per_file algorithm signal"}
folder=${3:?"min_size max_size path_instances nof_perms_per_file algorithm signal"}
ninst=${4:?"min_size max_size path_instances nof_perms_per_file algorithm signal"}
alg=${5:?"min_size max_size path_instances nof_perms_per_file algorithm signal"}
sig=${6:?"min_size max_size path_instances nof_perms_per_file algorithm signal"}


if [ ! -f "prog" ]; then
    cd ../../
    make
    make clean
    mv prog results/approx-big
    cd results/approx-big
fi

if [ ! -d "results-$alg" ]; then
    mkdir results-$alg
fi

for ((j=$minsize; j<=$maxsize; j+=5)); do
    fname="$folder/$j.txt"
    ./prog -a $alg -n $j -x 0 -s $sig -q $ninst -i $fname -v 1 >> $alg.out
    mv $alg.$j results-$alg
done

