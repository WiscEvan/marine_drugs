#!/usr/bin/env python

import pandas as pd


species_info = "/media/bigdrive2/evan/marine_drugs/marine_drugs/data/work/04/f13222453ffddc7ffbd2ad0cb7aaaf/busco_downloads/lineages/metazoa_odb10/info/species.info"
species_df = pd.read_csv(species_info,
    sep='\t', 
    header=None, 
    names=["OGID",'species_name'], 
    index_col='OGID')
#in species.info -> 400682	Amphimedon queenslandica
ogs_ids = "/media/bigdrive2/evan/marine_drugs/marine_drugs/data/work/04/f13222453ffddc7ffbd2ad0cb7aaaf/busco_downloads/lineages/metazoa_odb10/info/ogs.id.info"
dff = pd.read_csv(ogs_ids,
    sep='\t',
    header=None,
    names=["OGSID",'Busco id'])
# Each Busco ID has a HMM file of the form "<busco_id>.hmm"

dff["species_id"] = dff.OGSID.map(lambda x: x.split(":")[0])
dff["species"] = dff.species_id.map(lambda x: x.split("_")[0])

_species = df[df.species == "Amphimedon queenslandica"].iloc[0].name.astype(str)
dff[dff.species == _species]
# List of BUSCO IDs corresponding to Amphimedon queenslandica
_df = dff[dff.species == _species]
_buscos = _df["Busco id"].to_list()

# Now that we have BUSCO IDs we need to find them in the full_table.tsv
full_table = "/media/bigdrive2/evan/marine_drugs/marine_drugs/data/processed/sponge-markers/FL2015_4_markers/run_metazoa_odb10/full_table.tsv"

full_df = pd.read_csv(
    full_table, 
    sep='\t',
    header=None,
    names=["Busco id","Status","Sequence","Score","Length"], 
    comment="#")



# predictions.gff
# line == NODE_ ... AUGUSTUS ..
predictions = "/media/bigdrive2/evan/marine_drugs/marine_drugs/data/processed/sponge-markers/FL2015_4.predictions.gff"
contigs = []
with open(predictions) as fh:
    for line in fh:
        if "NODE" in line and "AUGUSTUS\ttranscript" in line:
            contig,*_,gene = line.strip().split("\t")
            contigs.append({"contig":contig, "gene":gene})

df = pd.DataFrame(contigs)
# predictions.aa
# >g24747.t1

# Now build contig->orf->Busco id full table from merge
master = pd.merge(full_df, df, how='left', left_on="Sequence", right_on="gene")
# Now we need to annotate whether Busco id is amphimedon or not...
master = pd.merge(master, dff, how='left', left_on="Busco id", right_on="Busco id")


master.species = master.species.astype(int)
master = pd.merge(master, species_df, how='left', left_on='species', right_on='OGID')

_species = df[df.species == "Amphimedon queenslandica"].iloc[0].name.astype(str)
master[master.species_name == "Amphimedon queenslandica"]
# List of BUSCO IDs corresponding to Amphimedon queenslandica

species_df = dff[dff.species == _species]