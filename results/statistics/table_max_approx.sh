###########################################################
# Created: Sex 28 Ago 2015 10:11:57 BRT
#
# Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
#
###########################################################
 



algs=("\\\texttt{2-PR}" "\\\texttt{2-PRx}" "\\\texttt{2-PSR}" "\\\texttt{2-PSRx}" "\\\texttt{2-PT}" "\\\texttt{2-PTx}" "\\\texttt{2-PST}" "\\\texttt{2-PSTx}" "\\\texttt{2-PRT}" "\\\texttt{2-PRTx}" "\\\texttt{2-PSRT}" "\\\texttt{2-PSRTx}" "\\\texttt{2-P{\RR}}" "\\\texttt{2-P{\RR}x}" "\\\texttt{2-PS{\RR}}" "\\\texttt{2-PS{\RR}x}" "\\\texttt{2-P{\RR}T}" "\\\texttt{2-P{\RR}Tx}" "\\\texttt{2-PS{\RR}T}" "\\\texttt{2-PS{\RR}Tx}")
files=("2-PR" "2-PRx" "2-PSR" "2-PSRx" "2-PT" "2-PTx" "2-PST" "2-PSTx" "2-PRT" "2-PRTx" "2-PSRT" "2-PSRTx" "2-SPR" "2-SPRx" "2-SPSR" "2-SPSRx" "2-SPRT" "2-SPRTx" "2-SPSRT" "2-SPSRTx")




echo "Small:"
for ((i=0; i<${#algs[@]}; i++)); do
    alg=${algs[$i]}
    file=${files[$i]}
    echo -e "${alg}\c"
    cat ../approx-small/${file}.approx | while read n avg max min; do 
        x=$(echo "scale = 3; $max/1" | bc)
        echo -e " & $x\c"
    done
    echo " \\\\ \\hline"
done

echo
echo "Big:"
for ((i=0; i<${#algs[@]}; i++)); do
    alg=${algs[$i]}
    file=${files[$i]}
    echo -e "${alg}\c"
    cat ../approx-big/${file}.approx | while read n avg max min; do 
        if [ "$n" == 25 ] || [ "$n" == 50 ] || [ "$n" == 75 ] || [ "$n" == 100 ] || [ "$n" == 125 ] || [ "$n" == 250 ] || [ "$n" == 500 ] || [ "$n" == 750 ] || [ "$n" == 1000 ] ; then
            x=$(echo "scale = 3; $max/1" | bc)
            echo -e " & $x\c"
        fi
    done
    echo " \\\\ \\hline"
done


numbers=(10 15 20 25 30 35 40 45 50 100 150 200 250 300 350 400 450 500 600 700 800 900 1000)
echo
echo "WPR:"
for ((a=2; a<=10; a++)); do
    echo "alpha = $a:"
    cat ../approx-small/WPR_$a.approx | while read n avg max min; do 
        x=$(echo "scale = 3; $max/1" | bc)
        echo -e " & $x\c"
    done

    cat ../approx-big/WPR_$a.approx | while read n avg max min; do 
        if [[ " ${numbers[@]} " =~ " $n " ]]; then
            x=$(echo "scale = 3; $max/1" | bc)
            echo -e " & $x\c"
        fi
    done
    echo " \\\\ \\hline"
    echo
done
