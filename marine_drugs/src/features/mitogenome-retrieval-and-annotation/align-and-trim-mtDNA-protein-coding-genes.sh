#!/usr/bin/env bash
#SBATCH --partition=queue
#SBATCH -N 1 # Nodes
#SBATCH -n 1 # Tasks
#SBATCH --cpus-per-task=1
#SBATCH --error=sponge_mtdna_phylogenetics.%J.err
#SBATCH --output=sponge_mtdna_phylogenetics.%J.out

DATA_DIR="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/data"
OUTDIR="${DATA_DIR}/processed/mitogenomes"

if [ ! -d $OUTDIR ];then
    mkdir -p $OUTDIR
fi

mitos_results_outdir="${DATA_DIR}/interim/mitogenomes/mitos_results"
mitos_parsed_results="${DATA_DIR}/interim/mitogenomes/mitos_results/munged"
# Retrieve mitos-annotated CDS seqs
# python /media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/src/features/mitogenome-retrieval-and-annotation/bin/parse_mitos_results.py \
#     --mitos $mitos_results_outdir \
#     --outdir $mitos_parsed_results

# 0. Get CDS seqs from lavrov
# NOTE: get_lavrov_et_al_mtdna.sh runs extract_seqs_from_gbk.py on downloaded gbk files
bash /media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/src/data/get_lavrov_et_al_mtdna.sh
# 0. Get CDS seqs from plese
# NOTE: get_plese_et_al_mtdna.sh runs translate_mtdna_cds_fna.py on downloaded mtdna.cds.fna files
bash /media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/src/data/get_plese_et_al_mtdna.sh


# 1.a Concatenate & sort individual CDS seqs
this_study_mtdna="${DATA_DIR}/interim/mitogenomes/mitos_results"
plese_mtdna="${DATA_DIR}/external/plese_et_al_2021"
lavrov_mtdna="${DATA_DIR}/external/lavrov_et_al"
concat_script="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/src/features/mitogenome-retrieval-and-annotation/bin/concat_all_mtdna_seqs.py"
concat_dir="${DATA_DIR}/interim/mitogenomes/concat" # concat {faa-dir}/{gene}.faa files
python $concat_script \
    --faa-dir $plese_mtdna $lavrov_mtdna $this_study_mtdna \
    --faa-format plese lavrov mitos \
    --out $concat_dir

# 1.b Check genes for minimum organism set
clean_dir="${DATA_DIR}/interim/mitogenomes/cleaned" # *.faa files
clean_script="/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/src/features/mitogenome-retrieval-and-annotation/bin/get_gene_superset.py"
python $clean_script \
    --in $clean_dir \
    --out $clean_dir

# 2. Perform (sorted) seq alignments using mafft
trim_dir="${DATA_DIR}/interim/mitogenomes/mafft" # *.mafft.faa
if [ ! -d $trim_dir ];then
    mkdir -p $trim_dir
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
    tree mytree = (Geodia_atlantica:0.0398506601,((((((((Cymbaxinella_damicornis:0.0427428194,(Stylissa_carteri:0.0000010039,Axinella_corrugata:0.0035119376):0.0025222437):0.0354941216,Agelas_schmidtii:0.1058205797):0.0400429958,(Halichondria_panicea:0.1482437273,(((((Phorbas_tenacior:0.0426456670,Phorbas_aerolatus:0.0590021411):0.0718548131,(FL2015_42:0.0000010039,FL2015_43:0.0000010039):0.1189112372):0.0251987687,Iotrochota_birotulata:0.0726808349):0.1178188342,(Cliona_varians:0.0419228182,((FL2015_5:0.0000000000,FL2015_4:0.0000000000):0.0000010039,FL2015_37:0.0000010039):0.0107642738):0.1070369602):0.0149193996,Tethya_actinia:0.1057451011):0.0119012796):0.0326320172):0.0114529970,(Phakellia_ventilabrum:0.0801861171,((Ptilocaulis_walpersi:0.0680416950,Ectyoplasia_ferox:0.0838579643):0.1415459585,Topsentia_ophiraphidites:0.0656570657):0.0193454332):0.0095672362):0.0068314882,((((((Dendrilla_antartica:0.0228991438,Igernella_notabilis:0.0417042420):0.3240517339,((((Spongia_officinalis:0.0095928627,(FL2014_9:0.0000010039,FL2014_3:0.0000010039):0.0088933133):0.0023914268,Hippospongia_lachne:0.0135319959):0.0065422855,Vaceletia_sp__GW948:0.0503482579):0.0035588816,(Ircinia_fasciculata:0.0075826306,(((FL2015_30:0.0000000000,FL2015_9:0.0000000000):0.0000010039,FL2015_34:0.0000010039):0.0004748177,FL2015_8:0.0004786944):0.0045377310):0.0155949492):1.1041240255):0.1900671837,(((Halisarca_caerulea:0.0214299731,Halisarca_dujardinii:0.0458689756):0.0124395984,Chondrilla_aff__nucula_CHOND:0.0229666497):0.1045185251,Aplysina_fulva:0.1238278208):0.0318670619):0.0252008185,((Plakinastrella_cf__onkodes_DVL_2011:0.2067712779,Oscarella_carmela:0.0965935083):0.0446049036,((Podila_verticillata:0.3974076051,Rhizopus_arrhizus:0.2617759872):0.2196543630,Allomyces_macrogynus:0.5523030356):0.7615065240):0.0771988612):0.0361889123,((Xestospongia_muta:0.0449223546,Callyspongia_plicifera:0.0850200249):0.0269989736,(Amphimedon_compressa:0.0013260560,FL2015_44:0.0004408276):0.1359561023):0.0166538627):0.0241673187,Ephydatia_muelleri:0.0722022123):0.0124062749):0.1323825196,Cinachyrella_kuekenthali:0.0504031530):0.0229253808,Poecillastra_compressa:0.0460104786):0.0091146516,Geodia_neptuni:0.0312706698):0.0059213994,Stryphnus_fortis:0.0642023557);
end;

begin mrbayes;

    [This block defines several different character sets that could be used in partitioning these data
    and then defines and enforces a partition called combined.
    NOTE: charset numbers were updated usign combined.nex]

    charset cob = 1-379;
    charset atp9 = 380-455;
    charset atp8 = 456-503;
    charset cox3 = 504-761;
    charset atp6 = 762-997;
    charset nad5 = 998-1566;
    charset nad2 = 1567-1950;
    charset cox2 = 1951-2183;
    charset nad3 = 2184-2299;
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
