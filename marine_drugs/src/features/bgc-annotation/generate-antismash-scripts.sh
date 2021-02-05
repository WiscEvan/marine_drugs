#!/usr/bin/env bash

SRC_DIR="$HOME/marine_drugs/marine_drugs/src/features/bgc-annotation"
template="${SRC_DIR}/antismash_template.sh"
sponges=(FL2014_9 FL2015_43 FL2015_44 FL2015_4 FL2015_5 FL2015_8 FL2015_9 FL2015_30 FL2015_34 FL2015_37 FL2015_42)
for sponge in ${sponges[@]};do
    script="${SRC_DIR}/antismash_${sponge}.sh"
    # Copy antismash_template.sh to antismash_FL2015_4.sh
    cp ${template} ${script}
    # edit antismash_FL2015_4.sh inputs and outputs filepaths
    sed -i "s,sample,${sponge},g" ${script}
    echo "wrote: ${script}"
done