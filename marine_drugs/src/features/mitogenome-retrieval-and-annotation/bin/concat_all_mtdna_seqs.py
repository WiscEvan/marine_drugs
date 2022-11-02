#!/usr/bin/env python

import argparse
import os
import re
import glob
from Bio import SeqIO
from typing import Literal, Dict,List
from functools import partial

from dataclasses import dataclass


@dataclass
class Gene:
    id: str
    description: str
    record: SeqIO.SeqRecord
    # organism: Optional[str]
    # accession: Optional[str]
    # genbank: Optional[str]

    def get_gene_name(self, format: Literal['mitos','lavrov','plese']) -> str:
        # Will first try to parse from lavrov format (shorthand methods)
        # Then will attempt to parse from mitos format
        # Finally plese format
        # Return None otherwise
        if format == "mitos":
            gene_name_getters = [self.get_mitos_gene_annotation]
        elif format == "lavrov":
            gene_name_getters = [
                self.get_nadh_shorthand,
                self.get_cob_shorthand,
                self.get_cox_shorthand,
                self.get_atp_shorthand,
                self.get_ribosomal_protein_shorthand,
                self.get_cantharalleus_literal,
            ]
        elif format == "plese":
            gene_name_getters = [self.get_plese_gene_annotation]
        else:
            raise ValueError(format)
        for gene_name_getter in gene_name_getters:
            gene_name = gene_name_getter()
            if gene_name:
                return gene_name.strip()

    def get_gene_shorthand(
        self, pattern: re.Pattern, group_num: int, gene: str = None
    ) -> str:
        match = pattern.search(self.description)
        if match:
            group = match.group(group_num)
            if gene:
                return f"{gene}{group}"
            else:
                return group

    def get_nadh_shorthand(self) -> str:
        # NADH dehydrogenase subunit : nad\d[L]
        pattern = re.compile(r"\S+\s-\sNADH dehydrogenase subunit\s(\S+)\s-\s(\S+)")
        gene_name = self.get_gene_shorthand(
            pattern=pattern,
            group_num=1,
            gene="nad",
        )
        if gene_name:
            return gene_name
        pattern = re.compile(r"\S+\s-\sNADH dehydrogenase, subunit\s(\S+)\s-\s(\S+)")
        return self.get_gene_shorthand(
            pattern=pattern,
            group_num=1,
            gene="nad",
        )

    def get_cob_shorthand(self) -> str:
        # apocytochrome b
        # cytochrome b
        pattern = re.compile(r"\S+\s-\sapocytochrome b\s-\s(\S+)")
        is_match = self.get_gene_shorthand(
            pattern=pattern,
            group_num=1,
        )
        if is_match:
            return "cob"
        pattern = re.compile(r"\S+\s-\scytochrome b\s-\s(\S+)")
        is_match = self.get_gene_shorthand(
            pattern=pattern,
            group_num=1,
        )
        if is_match:
            return "cob"

    def get_cox_shorthand(self) -> str:
        # cytochrome c oxidase subunit
        # cytochrome oxidase subunit 2
        # cytochrome c oxidase subunit III
        # cytochrome c oxidase subunit 3
        pattern = re.compile(r"\S+\s-\scytochrome c oxidase subunit\s(\S+)\s-\s(\S+)")
        raw_cyto_c_gene = self.get_gene_shorthand(
            pattern=pattern,
            group_num=1,
        )
        pattern = re.compile(r"\S+\s-\scytochrome oxidase subunit\s(\S+)\s-\s(\S+)")
        raw_cox_gene = self.get_gene_shorthand(
            pattern=pattern,
            group_num=1,
        )
        if not raw_cox_gene and not raw_cyto_c_gene:
            return

        if raw_cox_gene and "I" in raw_cox_gene:
            subunit_num = raw_cox_gene.count("I")
        elif raw_cyto_c_gene and "I" in raw_cyto_c_gene:
            subunit_num = raw_cyto_c_gene.count("I")
        elif raw_cox_gene:
            subunit_num = raw_cox_gene
        elif raw_cyto_c_gene:
            subunit_num = raw_cyto_c_gene
        return f"cox{subunit_num}"

    def get_atp_shorthand(self) -> str:
        # ATP synthase F0 subunit 8 : atp\d
        # atp9
        pattern = re.compile(r"\S+\s-\sATP synthase F0 subunit\s(\S+)\s-\s(\S+)")
        match = self.get_gene_shorthand(
            pattern=pattern,
            group_num=1,
            gene="atp",
        )
        if match:
            return match
        pattern = re.compile(r"\S+\s-\sATP synthetase subunit\s(\S+)\s-\s(\S+)")
        match = self.get_gene_shorthand(
            pattern=pattern,
            group_num=1,
            gene="atp",
        )
        if match:
            return match
        pattern = re.compile(r"\S+\s-\sH\(\+\)-transporting ATPase, F0 subunit\s(\S+)\s-\s(\S+)")
        return self.get_gene_shorthand(
            pattern=pattern,
            group_num=1,
            gene="atp",
        )

    def get_ribosomal_protein_shorthand(self) -> str:
        # ribosomal protein S8 : rpS8
        pattern = re.compile(r"\S+\s-\sribosomal protein\s(\S+)\s-\s(\S+)")
        return self.get_gene_shorthand(
            pattern=pattern,
            group_num=1,
            gene="rp",
        )
    
    def get_cantharalleus_literal(self) -> str:
        if "Cantharellus" in self.description:
            # Cantharellus formatted (is singleton from lavrov)
            return self.description.split('mt')[-1].split(';')[0].strip()

    def get_mitos_gene_annotation(self) -> str:
        if ";" in self.description:
            return self.description.split(';')[-1].strip().replace("_", "").replace("3'partial", "").replace('partial', "")
    
    def get_plese_gene_annotation(self) -> str:
        # >Haliclona_tubifera; 12044-13315; +; nad4_partial
        # Haliclona_tubifera; - 8119; +; nad5 3'partial - nad5 3'partial
        # Haliclona_tubifera; 9527-10333; +; cox3 3'partial
        # Haliclona_tubifera_atp6 3'partial
        # Haliclona_tubifera_cox2
        if "Code" in self.description:
            # Geodia_atlantica_nad3_Code - 4
            gene = self.description.split('_')[2]
        elif "Haliclona - tubifera" in self.description:
            # Haliclona - tubifera 9527-10333 + cox3 3'partial
            # Haliclona - tubifera 4210-4605 + nad1partial
            gene = self.description.split(' + ')[-1].replace("partial", "").replace("3'partial", "").strip()
        # >Haliclona_tubifera; 6206-8119; +; nad5 3'partial
        elif "Haliclona_tubifera" in self.description:
            # Haliclona_tubifera - atp6 3'partial
            gene = self.description.split('-')[-1].replace("3'partial", "").strip()
        elif "Phakellia_ventilabrum" in self.description and self.description.count("-") == 2:
            # Phakellia_ventilabrum - 4768-4971 + atp8
            gene = self.description.split(' + ')[-1].strip()
        elif "-" in self.description:
            gene = self.description.split('-')[-1].strip()
        gene = 'nad4' if gene == 'ad4' else gene # Phorbas_aerolatus - ad4
        return gene

    def find_organism(self, format: Literal['plese','lavrov']) -> str:
        # lavrov format: >NP_150603.1 - cytochrome c oxidase subunit II - Limulus polyphemus
        # 'ABW83951.1 - NADH dehydrogenase subunit 4 - Plakinastrella cf. onkodes DVL-2011'
        # Haliclona_tubifera; 8624-8860; +; atp8
        # Cantharellus cibarius mt cob ; 388 aa
        # plese format: >Cymbaxinella_damicornis_atp8
        # Cantharellus formatted (is singleton)
        if format == "plese":
            if "Haliclona - tubifera" in self.description:
                # Haliclona - tubifera 10387-11535 + cob
                organism = "Haliclona_tubifera"
            elif "Code" in self.description:
                # Geodia_atlantica_nad4l_Code - 4
                organism = "Geodia_atlantica"
            else:
                organism = self.description.split()[0]
        elif format == "lavrov":
            if "Cantharellus" in self.description:
                organism = self.description.split('mt')[0].strip()
            elif self.description.count("-") == 2:
                organism = self.description.split('-')[-1].strip()
            elif "-" in self.description and self.description.count("-") == 3:
                # NP_043745.1 - H(+)-transporting ATPase, F0 subunit 9 - Allomyces macrogynus
                # ABW83951.1 - NADH dehydrogenase subunit 4 - Plakinastrella cf. onkodes DVL-2011
                organism = self.description.split(' - ', 2)[-1].strip()
                # organism = organism.split('-')[-1] if '-' in organism else organism
                # organism = self.description.split(' - ')[-1]
            elif "_" in self.description:
                organism = '_'.join(self.description.split('_')[:2])
            else:
                organism = None
        else:
            raise ValueError(format)
        return organism


def get_genes(indir: str, format: Literal['mitos','lavrov','plese']) -> List[Gene]:
    genes = []
    for faa in glob.glob(os.path.join(indir, "*.faa")):
        for record in SeqIO.parse(faa, "fasta"):
            gene = Gene(
                id=record.id,
                description=record.description,
                record=record,
            )
            gene.genbank = os.path.basename(faa).replace(".faa", "")
            gene.organism = gene.find_organism(format=format)
            gene.gene_name = gene.get_gene_name(format=format)
            genes.append(gene)
    return genes


def get_mitos_results_protein_seqs(mitos: str) -> List[Gene]:
    prots_search_string = os.path.join(mitos, "**", "result.faa")
    prots = [fp for fp in glob.glob(prots_search_string, recursive=True)]
    mitos_genes = []
    for fp in prots:
        job_identifier = os.path.basename(os.path.dirname(os.path.dirname(fp)))
        sponge_id = job_identifier.replace("_code4_refseq89_metazoa", "")
        for record in SeqIO.parse(fp, "fasta"):
            # record.id is contig_id
            # record.description: NODE_703_length_19394_cov_1408.958694; 4902-5081; +; atp8
            # protein_name = record.description.split(';')[-1].strip()
            mitos_gene = Gene(
                id=record.id,
                description=record.description,
                record=record,
            )
            mitos_gene.organism = sponge_id
            mitos_gene.gene_name = mitos_gene.get_gene_name(format='mitos')
            mitos_genes.append(mitos_gene)
    return mitos_genes

def sort_genes(mtdna: List[Gene]) -> Dict[str,List[Gene]]:
    gene_records = {}
    for gene in mtdna:
        if not gene.gene_name:
            print(f"failed to recover gene name from {gene.description}")
        if gene.gene_name in gene_records:
            gene_records[gene.gene_name].append(gene)
        else:
            gene_records[gene.gene_name] = [gene]
    return gene_records

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--faa-dir", help="path to mtdna dir", nargs="+", required=True)
    parser.add_argument("--faa-format", help="path to mtdna dir", nargs="+", choices=["mitos", "plese", "lavrov"],required=True)
    parser.add_argument("--out", help="Output directory path", required=True)
    args = parser.parse_args()
    # Retrieve All CDS seqs from lavrov et. al. 2008
    # lavrov_dir = "/media/BRIANDATA3/kwan-bioinformatics-home-evan/evan/marine_drugs/marine_drugs/data/external/lavrov_et_al"
    gene_getter_dispatcher = {
        "mitos": get_mitos_results_protein_seqs,
        "plese": partial(get_genes, format="plese"),
        "lavrov": partial(get_genes, format="lavrov"),
    }
    genes = []
    for faa_dir,faa_format in zip(args.faa_dir, args.faa_format):
        gene_getter = gene_getter_dispatcher[faa_format]
        faa_dir_genes = gene_getter(faa_dir)
        genes.extend(faa_dir_genes)
        print(f"recovered {len(genes)} genes from {faa_dir}")
    sorted_genes = sort_genes(genes)

    # len(lavrov) == 62
    # len(plese) == 16
    # len(this study) == 12
    # total expected == 90
    # print(sorted_genes.keys())
    total_organisms = len(set([record.organism for records in sorted_genes.values() for record in records]))
    print(f"total organisms found: {total_organisms}")
    # print(set([record.organism for records in sorted_genes.values() for record in records]))
    for gene_name,records in sorted_genes.items():
        organisms = set()
        if not gene_name:
            continue
        n_organisms = len(set([r.organism for r in records]))
        fl_sponges_count = len([r.organism for r in records if "FL20" in r.organism])
        if fl_sponges_count != 12:
            outfname = f"{gene_name}.fasta" # Refer to mtdna w/o all sponges as *.fasta
            # To easily differentiate using *.faa instead of *.fasta
        else:
            print(f"{gene_name} contains all 12 sponges")
            outfname = f"{gene_name}.faa"
        outfpath = os.path.join(args.out, outfname)
        seqrecords = []
        for gene in records:
            record = gene.record
            record.description = ""
            record.name = gene.organism
            record.id = gene.organism
            if gene.organism in organisms:
                # print(f"Skipping duplicate for {gene_name} : {gene.organism}")
                continue
            seqrecords.append(record)
            organisms.add(gene.organism)
        # print(f'{gene_name} {len(records)}')
        n_written = SeqIO.write(seqrecords, outfpath, 'fasta')
        # print(f"Wrote {n_written} to {outfpath} ({n_organisms} organisms represented)")
    # Get counts of each organism in the separately written fasta files
    # grep -h ">" data/interim/mitogenomes/concat/*.faa | sort  | uniq -c

if __name__ == '__main__':
    main()