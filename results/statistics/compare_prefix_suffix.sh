###########################################################
# Created: Sex 28 Ago 2015 13:02:17 BRT
#
# Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
#
###########################################################
 

pathpref=${1?"*PATH_RESULTS_PREF* path_results_suff bin? sig? fileout fileoutperms"}
pathsuf=${2?"path_results_pref *PATH_RESULTS_SUFF* bin? sig? fileout fileoutperms"}
bin=${3?"path_results_pref path_results_suff *BIN?* sig? fileout fileoutperms"}
sig=${4?"path_results_pref path_results_suff BIN? *SIG?* fileout fileoutperms"}
fout=${5?"path_results_pref path_results_suff bin? sig? *FILEOUT* fileoutperms"}
fperms=${6?"path_results_pref path_results_suff bin? sig? fileout *FILEOUTPERMS*"}

gcc -c ../../permutations.c
gcc -c ../../breakpoints.c
gcc -c ../../util.c
gcc compare_prefix_suffix.c *.o -o compare_prefix_suffix -lm
rm *.o

./compare_prefix_suffix $pathpref $pathsuf $bin $sig $fout $fperms
