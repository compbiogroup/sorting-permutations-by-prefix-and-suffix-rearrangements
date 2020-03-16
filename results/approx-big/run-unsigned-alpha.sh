###########################################################
# Created: Qua 26 Ago 2015 16:49:57 BRT
#
# Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
#
###########################################################
 
algs=("WPRm" "WPRg" "WPR" "WPTg" "WPT" "WPRTg" "WPRT" "WPSRg" "WPSR" "WPSTg" "WPST" "WPSRTg" "WPSRT")

minsize=${1:?"min_size max_size path_instances nof_perms_per_file alpha"}
maxsize=${2:?"min_size max_size path_instances nof_perms_per_file alpha"}
folder=${3:?"min_size max_size path_instances nof_perms_per_file alpha"}
ninst=${4:?"min_size max_size path_instances nof_perms_per_file alpha"}
alpha=${5:?"min_size max_size path_instances nof_perms_per_file alpha"}

if [ ! -f "prog" ]; then
    cd ../../
    make
    make clean
    mv prog results/approx-big
    cd results/approx-big
fi

for a in ${algs[@]}; do
    if [ ! -d "results-${a}_${alpha}" ]; then
        mkdir results-${a}_${alpha}
    fi
done

for ((j=$minsize; j<=$maxsize; j+=5)); do
    for a in ${algs[@]}; do
        fname="$folder/$j.txt"
        ./prog -a $a -n $j -x $alpha -q $ninst -s 0 -i $fname -v 1 >> ${a}_${alpha}.out
        mv ${a}*.$j results-${a}_${alpha}
    done
done
