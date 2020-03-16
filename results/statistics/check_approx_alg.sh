###########################################################
# Created: Sex 28 Ago 2015 11:08:02 BRT
#
# Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
#
###########################################################
 



pathresults=${1?"*PATH_RESULTS* alg_name minsize maxsize step fileout factor"}
alg=${2?"path_results *ALG_NAME* minsize maxsize step *FILEOUT* factor"}
minsize=${3?"path_results alg_name *MINSIZE* maxsize step fileout factor"}
maxsize=${4?"path_results alg_name minsize *MAXSIZE* step fileout factor"}
step=${5?"path_results alg_name minsize maxsize *STEP* fileout factor"}
fout=${6?"path_results alg_name minsize maxsize step *FILEOUT* factor"}
f=${7?"path_results alg_name minsize maxsize step fileout *FACTOR*"}

gcc check_approx_alg.c -o check_approx_alg

for ((n=$minsize; n<=$maxsize; n+=$step)); do
    ./check_approx_alg $n $pathresults/$alg.$n.approx $f >> $fout
done
