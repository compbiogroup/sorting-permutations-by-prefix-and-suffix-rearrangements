###########################################################
# Created: Sex 28 Ago 2015 10:11:57 BRT
#
# Author: Carla N. Lintzmayer, carlanl@ic.unicamp.br
#
###########################################################
 



algs=("2-PR" "2-PRx" "2-PSR" "2-PSRx" "2-PT" "2-PTx" "2-PST" "2-PSTx" "2-PRT" "2-PRTx" "2-PSRT" "2-PSRTx" "2-SPR" "2-SPRx" "2-SPSR" "2-SPSRx" "2-SPRT" "2-SPRTx" "2-SPSRT" "2-SPSRTx")
probs=("pr" "pr" "psr" "psr" "pt" "pt" "pst" "pst" "prt" "prt" "psrt" "psrt" "spr" "spr" "spsr" "spsr" "sprt" "sprt" "spsrt" "spsrt")
dists=("pr" "pr" "psr" "psr" "pt" "pt" "pst" "pst" "prt" "prt" "psrt" "psrt" "spr" "spr" "spsr" "spsr" "sprt" "sprt" "spsrt" "spsrt")
sigs=("0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "1" "1" "1" "1" "1" "1" "1" "1")
maxn=("10" "10" "10" "10" "10" "10" "10" "10" "10" "10" "10" "10" "9" "9" "9" "9" "9" "9" "9" "9")

for ((i=0; i<${#probs[@]}; i++)); do
    prob=${probs[$i]}
    alg=${algs[$i]}
    sig=${sigs[$i]}
    ./generate_approx.sh ../approx-big/results-$alg $alg 10 1000 5 T $prob 0 null 0 $sig ../approx-big/$alg.approx $alg
done


for ((i=0; i<${#probs[@]}; i++)); do
    prob=${probs[$i]}
    alg=${algs[$i]}
    dist=${dists[$i]}
    sig=${sigs[$i]}
    n=${maxn[$i]}
    ./generate_approx.sh ../approx-small/results-$alg $alg 2 $n 1 D $prob 0 ~/Workspace/distances/$dist 1 $sig ../approx-small/$alg.approx $alg
done

#
algs=("WPRm" "WPR" "WPRg" "WPSR" "WPSRg" "WPT" "WPTg" "WPST" "WPSTg" "WPRT" "WPRTg" "WPSRT" "WPSRTg" "WSPR" "WSPRg" "WSPSR" "WSPSRg" "WSPRT" "WSPRTg" "WSPSRT" "WSPSRTg")
probs=("wpr" "wpr" "wpr" "wpsr" "wpsr" "wpt" "wpt" "wpst" "wpst" "wprt" "wprt" "wpsrt" "wpsrt" "wspr" "wspr" "wspsr" "wspsr" "wsprt" "wsprt" "wspsrt" "wspsrt")
dists=("wpr" "wpr" "wpr" "wpsr" "wpsr" "wpt" "wpt" "wpst" "wpst" "wprt" "wprt" "wpsrt" "wpsrt" "wspr" "wspr" "wspsr" "wspsr" "wsprt" "wsprt" "wspsrt" "wspsrt")
sigs=("0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "0" "1" "1" "1" "1" "1" "1" "1" "1")
maxn=("10" "10" "10" "10" "10" "10" "10" "10" "10" "10" "10" "10" "10" "8" "8" "8" "8" "8" "8" "8" "8")

for ((i=0; i<${#probs[@]}; i++)); do
    prob=${probs[$i]}
    alg=${algs[$i]}
    sig=${sigs[$i]}
    ./generate_approx.sh ../approx-big/results-${alg}_1 ${alg}_1.000 10 1000 5 T $prob 1 null 0 $sig ../approx-big/${alg}_1.approx ${alg}_1.000
done

for ((i=0; i<${#probs[@]}; i++)); do
    prob=${probs[$i]}
    alg=${algs[$i]}
    dist=${dists[$i]}
    sig=${sigs[$i]}
    n=${maxn[$i]}
    ./generate_approx.sh ../approx-small/results-${alg}_1 ${alg}_1.000 2 $n 1 D $prob 1 ~/Workspace/distances/$dist 0 $sig ../approx-small/${alg}_1.approx ${alg}_1.000
done

#
for ((a=2; a<=10; a++)); do
    alg="WPR_${a}.000"
    prob="wpr"
    dist="wpr"
    sig=0
    n=10
    ./generate_approx.sh ../approx-small/results-WPR_${a} WPR_${a}.000 2 $n 1 D $prob $a ~/Workspace/distances/$dist.$a 0 $sig ../approx-small/WPR_${a}.approx WPR_${a}.000
done

for ((a=2; a<=10; a++)); do
    alg="WPR_${a}.000"
    prob="wpr"
    dist="wpr"
    sig=0
    ./generate_approx.sh ../approx-big/results-WPR_${a} WPR_${a}.000 10 1000 5 T $prob $a null 0 $sig ../approx-big/WPR_${a}.approx WPR_${a}.000
done
