#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=sponge_mtdna_phylogenetics.%J.err
#SBATCH --output=sponge_mtdna_phylogenetics.%J.out

DATA_DIR="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/data"
OUTDIR="${DATA_DIR}/processed/mitogenomes"

mitos_results_outdir="${DATA_DIR}/interim/mitogenomes/mitos_results"
mitos_parsed_results="${DATA_DIR}/interim/mitogenomes/mitos_results/munged"
# Retrieve mitos-annotated CDS seqs
# python /media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/src/features/mitogenome-retrieval-and-annotation/bin/parse_mitos_results.py \
#     --mitos $mitos_results_outdir \
#     --outdir $mitos_parsed_results

# 0. Get CDS seqs from external references (lavrov and plese) and mitos
bash /media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/src/data/get_lavrov_et_al_mtdna.sh
bash /media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/src/data/get_plese_et_al_mtdna.sh
bash /media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/src/data/extract_reference_mtdna_gbk_to_faa.sh

# 1.a Concatenate & sort individual CDS seqs
this_study_mtdna="${DATA_DIR}/interim/mitogenomes/mitos_results"
plese_mtdna="${DATA_DIR}/external/plese_et_al_2021"
lavrov_mtdna="${DATA_DIR}/external/lavrov_et_al"
concat_script="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/src/features/mitogenome-retrieval-and-annotation/bin/concat_all_mtdna_seqs.py"
outdir="${DATA_DIR}/interim/mitogenomes/concat"
python $concat_script \
    --faa-dir $plese_mtdna $lavrov_mtdna $this_study_mtdna \
    --faa-format plese lavrov mitos \
    --out $outdir

# 1.b Check genes for minimum organism set
concat_dir="${DATA_DIR}/interim/mitogenomes/concat" # *.faa files
clean_dir="${DATA_DIR}/interim/mitogenomes/cleaned" # *.faa files
clean_script="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/src/features/mitogenome-retrieval-and-annotation/bin/get_gene_superset.py"
python $clean_script \
    --in $concat_dir \
    --out $clean_dir

# 2. Perform (sorted) seq alignments using mafft
clean_dir="${DATA_DIR}/interim/mitogenomes/cleaned" # *.faa files
align_dir="${DATA_DIR}/interim/mitogenomes/mafft" # *.mafft.faa
if [ ! -d $align_dir ];then
    mkdir -p $align_dir
fi

# mamba install -c bioconda mafft
for faa in `find $clean_dir -name "*.faa"`;do
    echo "aligning $(basename ${faa})"
    outfname=$(basename ${faa/.faa/.mafft.faa})
    outfpath="${align_dir}/${outfname}"
    mafft --auto $faa > $outfpath
done

# 3. Trim alignments using BMGE
trim_dir="${DATA_DIR}/interim/mitogenomes/bmge" # *.bmge.nex
if [ ! -d $trim_dir ];then
    mkdir -p $trim_dir
fi
# mamba install -c bioconda bmge
for faa in `find $align_dir -name "*.mafft.faa"`;do
    echo "trimming $(basename ${faa})"
    outfname=$(basename ${faa/.mafft.faa/.bmge.nex})
    outfpath="${trim_dir}/${outfname}"
    bmge -i $faa -t AA -on $outfpath
done

# 4. Concatenate sorted, aligned and trimmed files (these are in nexus format output from bmge)
combine_nexi_script="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/src/features/mitogenome-retrieval-and-annotation/bin/concat_nex.py"
combined_nexus="${OUTDIR}/combined.nex"
if [ ! -d $OUTDIR ];then
    mkdir -p $OUTDIR
fi
nexi=(`find $trim_dir -name "*.nex"`)
python $combine_nexi_script \
    --nex ${nexi[@]} \
    --out $combined_nexus

# 5. Create maximum-likelihood tree with IQ-TREE
# mamba install -c bioconda iqtree -y
iqtree -nt AUTO -s $combined_nexus

# NOTE: This is an optional step that renames the tree terminals from FL201{4,5}_* to identified species
iqtree_treefile="${OUTDIR}/combined.nex.treefile"
rename_tree_terminals_script="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/src/features/mitogenome-retrieval-and-annotation/bin/rename_tree_sponge_ids.py"
python $rename_tree_terminals_script \
    --tree $iqtree_treefile \
    --out ${iqtree_treefile/.treefile/.renamed.treefile} \
    --metadata "${DATA_DIR}/raw/sponge_metadata.tsv"

# 6. Create bayesian inference tree using Mr.Bayes

# 6.1 Copy/Paste mrbayes_txt chunk below (with iqtree information) in to concatenated nexus file
# NOTE: Retrieve charset definitions from combined.nex
charset=$(grep -h -e "charset" $combined_nexus)
charpartition=$(grep -h -e "charpartition" $combined_nexus | sed -e 's/charpartition/partition/g')
# NOTE: 
#       - partition should follow mrbayes format, e.g.
#           partition combined = 9: cob, atp9, atp8, cox3, atp6, nad5, nad2, cox2, nad3;
#       - Copy/Paste $mrbayes_txt at the end of the $combined_nexus file prior to running mrbayes
#       - copy/paste contents of $iqtree_treefile to TREE block, e.g.
#           tree mytree = {paste tree here}
mrbayes_txt="""
begin trees;
    [NOTE: tree was copy/pasted from iqtree output]
    tree mytree = (Geodia_atlantica:0.0394863989,((((((((Cymbaxinella_damicornis:0.0427254052,(Stylissa_carteri:0.0000010023,Axinella_corrugata:0.0034820358):0.0024962957):0.0350022973,Agelas_schmidtii:0.1057058027):0.0403760317,(Halichondria_panicea:0.1477771616,(((((Phorbas_tenacior:0.0425429588,Phorbas_aerolatus:0.0584839434):0.0713840062,(FL2015_42:0.0000010023,FL2015_43:0.0000010023):0.1188023226):0.0254036761,Iotrochota_birotulata:0.0721890769):0.1177001660,(Cliona_varians:0.0418436541,((FL2015_5:0.0000000000,FL2015_4:0.0000000000):0.0000010023,FL2015_37:0.0000010023):0.0107409308):0.1069456874):0.0147042238,Tethya_actinia:0.1054835227):0.0118103667):0.0333500062):0.0118973080,(Phakellia_ventilabrum:0.0802092926,((Ptilocaulis_walpersi:0.0682455487,Ectyoplasia_ferox:0.0841239407):0.1401877706,Topsentia_ophiraphidites:0.0658592237):0.0193997120):0.0095500239):0.0067864669,((((((Dendrilla_antartica:0.0229190468,Igernella_notabilis:0.0414288423):0.3233820559,((((Spongia_officinalis:0.0092685695,(FL2014_9:0.0000010023,FL2014_3:0.0000010023):0.0092401240):0.0017211796,Hippospongia_lachne:0.0142005040):0.0063368838,(((FL2015_30:0.0000000000,FL2015_9:0.0000000000):0.0000010023,FL2015_34:0.0000010023):0.0000010023,FL2015_8:0.0009535481):0.0233645928):0.0054098786,Vaceletia_sp__GW948:0.0454024732):1.1023409658):0.1908932719,(((Halisarca_caerulea:0.0214182581,Halisarca_dujardinii:0.0451653834):0.0124480773,Chondrilla_aff__nucula_CHOND:0.0228502544):0.1033014045,Aplysina_fulva:0.1240344404):0.0319623781):0.0250786400,((Plakinastrella_cf__onkodes_DVL_2011:0.2066553555,Oscarella_carmela:0.0966001045):0.0445247603,((Podila_verticillata:0.3972800455,Rhizopus_arrhizus:0.2615529299):0.2201387721,Allomyces_macrogynus:0.5517200469):0.7649240317):0.0761084591):0.0371447353,((Xestospongia_muta:0.0438892667,Callyspongia_plicifera:0.0849954450):0.0267689885,(Amphimedon_compressa:0.0013327601,FL2015_44:0.0004433698):0.1360057297):0.0160186420):0.0237396436,Ephydatia_muelleri:0.0720468714):0.0124111438):0.1326045873,Cinachyrella_kuekenthali:0.0495809450):0.0227277527,Poecillastra_compressa:0.0457719397):0.0091010630,Geodia_neptuni:0.0314393112):0.0059688302,Stryphnus_fortis:0.0645336638);
end;

begin mrbayes;

    [This block defines several different character sets that could be used in partitioning these data
    and then defines and enforces a partition called combined.]

    charset cob = 1-379;
    charset atp9 = 380-455;
    charset atp8 = 456-504;
    charset cox3 = 505-762;
    charset atp6 = 763-998;
    charset nad5 = 999-1564;
    charset nad2 = 1565-1945;
    charset cox2 = 1946-2178;
    charset nad3 = 2179-2294;
    partition combined = 9: cob, atp9, atp8, cox3, atp6, nad5, nad2, cox2, nad3;
    set partition=combined;
    lset applyto=(all) nst=6 rates=invgamma;
    startvals tau=mytree;
    mcmcp ngen=15000000 relburnin=yes burninfrac=0.25 nchains=4 samplefreq=1000 Nruns=2 nperts=4 stoprule=yes stopval=0.01;
    [
        The following lines set up a particular model (the one discussed in the MrBayes manual). Uncomment this block
        if you want to set up the model when executing this file instead of specifying it yourself.
        i.e. execute /path/to/combined.nex
        # NOTE: unlink command resulted in error
        unlink statefreq=(all) revmat=(all) shape=(all) pinvar=(all);
        # NOTE:
        # run command: mcmc to start inference
        # After `mcmc` is finished (and has converged) run sumt and sump
        sumt burnin=250;
        sump burnin=250;
    ]

end;
"""

# 6.2 download & install mr bayes
# mb="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/src/features/mitogenome-retrieval-and-annotation/MrBayes/bin/mb"
# 6.3 Start mr bayes using command:
#   mb
# 6.4 Read in nexus file (with tree data and parameters)
#   execute </path/to/nexus/file>
# e.g. `Execute data/processed/mitogenomes/combined.nex`
# 6.5 Run Markov chain
#   mcmc
# ... will start run and stop at average std.dev. of split freq. (ASDSF) less than or equal to 0.01
