###########################################################
# Created: Mon Dec 17 15:38:05 2012
#
# Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
#
###########################################################

algs=("2-PSR" "2-PST" "2-PSRT")
sizes=(10 50 100 150 200 250)

mult=${1:?"mult path_instances"}
folder=${2:?"mult path_instances"}

if [ ! -f "prog" ]; then
    cd ../../
    make
    make clean
    mv prog results/approx-huge
    cd results/approx-huge
fi

for a in ${algs[@]}; do
    if [ ! -d "results-${a}_${mult}" ]; then
        mkdir results-${a}_${mult}
    fi
done

for n in ${sizes[@]}; do
    for a in ${algs[@]}; do
        fname="$folder/${n}.$((n*mult)).txt"
        ./prog -a $a -n $n -x 0 -s 0 -q $((n*mult)) -i $fname -v 1 >> ${a}_${mult}.out
        mv $a.$n results-${a}_${mult}
    done
done

