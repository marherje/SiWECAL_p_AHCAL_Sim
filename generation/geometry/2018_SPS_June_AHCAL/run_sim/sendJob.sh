#!/bin/zsh

nqueue=1
neventsPerProcess=100

#for particle (e- pi-)
for particle (mu-)
do
  #for energy (100 150)
  for energy (10)
  do
    cp run.sh run_${particle}_${energy}GeV.sh
    /bin/vi run_${particle}_${energy}GeV.sh <<EOF >& /dev/null
:%s/NEVENTSPERPROCESS/$neventsPerProcess/g
:w!
:q
EOF

    cp htc_jobSettings.sub htc_jobSettings_${particle}_${energy}GeV.sub
    /bin/vi htc_jobSettings_${particle}_${energy}GeV.sub <<EOF >& /dev/null
:%s/PARTICLE/${particle}/g
:%s/ENERGY/$energy/g
:%s/NQUEUE/$nqueue/g
:w!
:q
EOF
    condor_submit htc_jobSettings_${particle}_${energy}GeV.sub

  done
done
