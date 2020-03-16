###########################################################
# Created: Sex 28 Ago 2015 11:19:11 BRT
#
# Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
#
###########################################################
 

pathresults1=${1?"*PATH_RESULTS1* alg_name1 path_results2 alg_name2 minsize maxsize step fileout"}
alg1=${2?"path_results1 *ALG_NAME1* path_results2 alg_name2 minsize maxsize step *FILEOUT*"}
pathresults2=${3?"path_results1 alg_name1 *PATH_RESULTS2* alg_name2 minsize maxsize step fileout"}
alg2=${4?"path_results1 alg_name1 path_results2 *ALG_NAME2* minsize maxsize step *FILEOUT*"}
minsize=${5?"path_results1 alg_name1 path_results2 alg_name2 *MINSIZE* maxsize step fileout"}
maxsize=${6?"path_results1 alg_name1 path_results2 alg_name2 minsize *MAXSIZE* step fileout"}
step=${7?"path_results1 alg_name1 path_results2 alg_name2 minsize maxsize *STEP* fileout"}
fout=${8?"path_results1 alg_name1 path_results2 alg_name minsize maxsize step *FILEOUT*"}

gcc compare_two_algs.c -o compare_two_algs

for ((n=$minsize; n<=$maxsize; n+=$step)); do
    ./compare_two_algs $n $pathresults1/$alg1.$n $pathresults2/$alg2.$n >> $fout
done
