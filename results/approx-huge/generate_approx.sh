###########################################################
# Created: Qua 26 Ago 2015 16:01:01 BRT
#
# Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
#
###########################################################
 


pathresults=${1?"path_results alg_file_name problem_name sig file_out_total_name file_out_perms_name"}
alg=${2?"path_results alg_file_name problem_name sig file_out_total_name file_out_perms_name"}
problem=${3?"path_results alg_file_name problem_name sig file_out_total_name file_out_perms_name"}
sig=${4?"path_results alg_file_name problem_name sig file_out_total_name file_out_perms_name"}
fout=${5?"path_results alg_file_name problem_name sig file_out_total_name file_out_perms_name"}
foutname=${6?"path_results alg_file_name problem_name sig file_out_total_name file_out_perms_name"}

sizes=(10 50 100 150 200 250)

gcc -c ../../breakpoints.c
gcc -c ../../permutations.c
gcc -c ../../util.c
gcc -c ../../rearrangements.c
gcc generate_approx.c *.o -o generate_approx -lm
rm *.o

for n in ${sizes[@]}; do
    fresults=$pathresults/$alg.$n
    foutind=$pathresults/$foutname.$n.approx
    ./generate_approx $fresults $n $problem $sig $foutind >> $fout
done

