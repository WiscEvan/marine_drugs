#!/usr/bin/env bash


for i in chocophlan,full uniref,uniref50_diamond utility_mapping,full;do
    IFS=",";
    set -- $i;
    echo "downloading database: $1 and $2 (via Docker)";
    docker run \
        -u $(id -g):$(id -u) \
        --rm \
        --detach=TRUE \
        -v /media/bigdrive1/Databases/humann_databases/:/databases:rw \
        biobakery/humann:latest \
            humann_databases --download $1 $2 /databases/
done