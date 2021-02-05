
## AntiSMASH annotation for sponge metagenomes

Naive running with AntiSMASH wrapper and docker container

```bash
bash generate-antismash-scripts.sh
sponges=(FL2014_3 FL2014_9 FL2015_43 FL2015_44 FL2015_4 FL2015_5 FL2015_8 FL2015_9 FL2015_30 FL2015_34 FL2015_37 FL2015_42)
for sponge in ${sponges[@]};do
    sbatch antismash_${songe}.sh
done
```

Could specify the docker container and run samples via nextflow. This may be more fault tolerant?