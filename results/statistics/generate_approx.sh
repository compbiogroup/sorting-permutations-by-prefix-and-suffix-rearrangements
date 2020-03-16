###########################################################
# Created: Qua 26 Ago 2015 16:01:01 BRT
#
# Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
#
###########################################################
 


pathresults=${1?"*PATH_RESULTS* alg_file_name minsize maxsize step type_lower_bound(D|T) problem_name_for_lower_bound alpha path_dists bin sig file_out_total_name file_out_perms_name"}
alg=${2?"path_results *ALG_FILE_NAME* minsize maxsize step type_lower_bound(D|T) problem_name_for_lower_bound alpha path_dists bin sig file_out_total_name file_out_perms_name"}
minsize=${3?"path_results alg_file_name *MINSIZE* maxsize step type_lower_bound(D|T) problem_name_for_lower_bound alpha path_dists bin sig file_out_total_name file_out_perms_name"}
maxsize=${4?"path_results alg_file_name minsize *MAXSIZE* step type_lower_bound(D|T) problem_name_for_lower_bound alpha path_dists bin sig file_out_total_name file_out_perms_name"}
step=${5?"path_results alg_file_name minsize maxsize *STEP* type_lower_bound(D|T) problem_name_for_lower_bound alpha path_dists bin sig file_out_total_name file_out_perms_name"}
lb=${6?"path_results alg_file_name minsize maxsize step *TYPE_LOWER_BOUND(D|T)* problem_name_for_lower_bound alpha path_dists bin sig file_out_total_name file_out_perms_name"}
problem=${7?"path_results alg_file_name minsize maxsize step type_lower_bound(D|T) *PROBLEM_NAME_FOR_LOWER_BOUND* alpha path_dists bin sig file_out_total_name file_out_perms_name"}
alpha=${8?"path_results alg_file_name minsize maxsize step type_lower_bound(D|T) problem_name_for_lower_bound *ALPHA* path_dists bin sig file_out_total_name file_out_perms_name"}
pathdists=${9?"path_results alg_file_name minsize maxsize step type_lower_bound(D|T) problem_name_for_lower_bound alpha *PATH_DISTS* bin sig file_out_total_name file_out_perms_name"}
bin=${10?"path_results alg_file_name minsize maxsize step type_lower_bound(D|T) problem_name_for_lower_bound alpha path_dists *BIN* sig file_out_total_name file_out_perms_name"}
sig=${11?"path_results alg_file_name minsize maxsize step type_lower_bound(D|T) problem_name_for_lower_bound alpha path_dists bin *SIG* file_out_total_name file_out_perms_name"}
fout=${12?"path_results alg_file_name minsize maxsize step type_lower_bound(D|T) problem_name_for_lower_bound alpha path_dists bin sig *FILE_OUT_TOTAL_NAME* file_out_perms_name"}
foutname=${13?"path_results alg_file_name minsize maxsize step type_lower_bound(D|T) problem_name_for_lower_bound alpha path_dists bin sig file_out_total_name *FILE_OUT_PERMS_NAME"}

gcc -c ../../breakpoints.c
gcc -c ../../permutations.c
gcc -c ../../util.c
gcc -c ../../rearrangements.c
gcc generate_approx.c *.o -o generate_approx -lm
rm *.o

for ((n=$minsize; n<=$maxsize; n+=$step)); do
    fresults=$pathresults/$alg.$n
    fdists=$pathdists/$n.txt
    foutind=$pathresults/$foutname.$n.approx
    ./generate_approx $fresults $n $lb $problem $alpha $fdists $bin $sig $foutind >> $fout
done

