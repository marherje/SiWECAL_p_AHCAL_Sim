mkdir LCFIPLOTS
mkdir LOGS
for energy in 250 500
do
    mkdir LCFIPLOTS/${energy} 
    for pol in eL_pR eR_pL
    do
	mkdir LCFIPLOTS/${energy}/${pol}
	for cat in A B C D
	do
	    mkdir LCFIPLOTS/${energy}/${pol}/${cat}
	    root -l -q lcfiplusvariables.cc\(\"${energy}\",\"${pol}\",\"${cat}\"\) > log_${energy}${pol}_${cat}
	    mv *png LCFIPLOTS/${energy}/${pol}/${cat}/.
	    mv log_${energy}${pol}_${cat} LOGS/.
	done	
    done
done
