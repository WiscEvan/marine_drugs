#!/usr/bin/env bash

outdir="$HOME/sponge_paper/sponge_paper/data/raw/assemblies"

declare -a samples
samples=(FL2014_3 FL2014_9 FL2015_42 FL2015_43 FL2015_9 FL2015_30 FL2015_34 FL2015_4 FL2015_5 FL2015_37 FL2015_8)

for sample in ${samples[@]};do
    tarball="/media/external2/Sponge_assemblies/${sample}_Spades_output.tar.gz"
    listing="${outdir}/${sample}_spades_tarball_listing.txt"
    if [ ! -f $listing ];then
        echo "listing ${tarball} contents to ${listing}:"
        tar -tf ${tarball} > ${listing}
    fi
done

# tar -tf /media/external2/Sponge_assemblies/FL2014_3_Spades_output.tar.gz > "${outdir}/FL2014_3_Spades_output_tarball_listing.txt"
# tar -tf /media/external2/Sponge_assemblies/FL2014_9_Spades_output.tar.gz > "${outdir}/FL2014_9_Spades_output_tarball_listing.txt"
# tar -tf /media/external2/Sponge_assemblies/FL2015_42_Spades_output.tar.gz > "${outdir}/FL2015_42_Spades_output_tarball_listing.txt"
# tar -tf /media/external2/Sponge_assemblies/FL2015_43_Spades_output.tar.gz > "${outdir}/FL2015_43_Spades_output_tarball_listing.txt"
# tar -tf /media/external2/Sponge_assemblies/FL2015_9_Spades_output.tar.gz > "${outdir}/FL2015_9_Spades_output_tarball_listing.txt"
# tar -tf /media/external2/Sponge_assemblies/FL2015_30_Spades_output.tar.gz > "${outdir}/FL2015_30_Spades_output_tarball_listing.txt"
# tar -tf /media/external2/Sponge_assemblies/FL2015_34_Spades_output.tar.gz > "${outdir}/FL2015_34_Spades_output_tarball_listing.txt"
# tar -tf /media/external2/Sponge_assemblies/FL2015_4_Spades_output.tar.gz > "${outdir}/FL2015_4_Spades_output_tarball_listing.txt"
# tar -tf /media/external2/Sponge_assemblies/FL2015_5_Spades_output.tar.gz > "${outdir}/FL2015_5_Spades_output_tarball_listing.txt"
# # tar -tf /media/external2/Sponge_assemblies/FL2015_44_Spades_output.tar.gz > "${outdir}/FL2015_44_Spades_output_tarball_listing.txt" This is a failed assembly
# tar -tf /media/external2/Sponge_assemblies/FL2015_37_Spades_output.tar.gz > "${outdir}/FL2015_37_Spades_output_tarball_listing.txt"
# tar -tf /media/external2/Sponge_assemblies/FL2015_8_Spades_output.tar.gz > "${outdir}/FL2015_8_Spades_output_tarball_listing.txt"