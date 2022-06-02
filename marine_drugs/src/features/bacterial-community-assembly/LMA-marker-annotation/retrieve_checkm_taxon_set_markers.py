#!/usr/bin/env python

import argparse
import os
import subprocess


def read_checkm_taxon_set(line):
    """Construct bin marker set data from line.

    NOTE
    ----

        Class and func yoinked from:
        https://github.com/Ecogenomics/CheckM/blob/e79d420194f385c67d3b3710c6beadcdf710598a/checkm/markerSets.py#L157
    
    """
    markerSets = []
    lineSplit = line.split("\t")
    numMarkerSets = int(lineSplit[1])
    for i in range(0, numMarkerSets):
        uid = lineSplit[i * 4 + 2]
        lineageStr = lineSplit[i * 4 + 3]
        numGenomes = int(lineSplit[i * 4 + 4])
        markerSet = eval(lineSplit[i * 4 + 5])
        markerSets.append(MarkerSet(uid, lineageStr, numGenomes, markerSet))
    return markerSets


class MarkerSet:
    """A collection of marker genes organized into co-located sets."""

    def __init__(self, UID, lineageStr, numGenomes, markerSet):

        self.UID = UID  # unique ID of marker set
        self.lineageStr = lineageStr  # taxonomic string associated with marker set
        self.numGenomes = numGenomes  # number of genomes used to calculate marker set
        self.markerSet = markerSet  # marker genes organized into co-located sets

    def __repr__(self):
        return (
            str(self.UID)
            + "\t"
            + self.lineageStr
            + "\t"
            + str(self.numGenomes)
            + "\t"
            + str(self.markerSet)
        )

    def size(self):
        """Number of marker genes and marker gene sets."""
        numMarkerGenes = 0
        for m in self.markerSet:
            numMarkerGenes += len(m)

        return numMarkerGenes, len(self.markerSet)

    def numMarkers(self):
        """Number of marker genes."""
        return self.size()[0]

    def numSets(self):
        """Number of marker sets."""
        return len(self.markerSet)

    def getMarkerGenes(self):
        """Get marker genes within marker set."""
        markerGenes = set()
        for m in self.markerSet:
            for marker in m:
                markerGenes.add(marker)

        return markerGenes

    def removeMarkers(self, markersToRemove):
        """Remove specified markers from marker sets."""
        newMarkerSet = []
        for ms in self.markerSet:
            newMS = ms - markersToRemove

            if len(newMS) != 0:
                newMarkerSet.append(newMS)

        self.markerSet = newMarkerSet

    def genomeCheck(self, hits, bIndividualMarkers):
        """Calculate genome completeness and contamination."""
        if bIndividualMarkers:
            present = 0
            multiCopyCount = 0
            for marker in self.getMarkerGenes():
                if marker in hits:
                    present += 1
                    multiCopyCount += len(hits[marker]) - 1

            percComp = 100 * float(present) / self.numMarkers()
            percCont = 100 * float(multiCopyCount) / self.numMarkers()
        else:
            comp = 0.0
            cont = 0.0
            for ms in self.markerSet:
                present = 0
                multiCopy = 0
                for marker in ms:
                    count = len(hits.get(marker, []))
                    if count == 1:
                        present += 1
                    elif count > 1:
                        present += 1
                        multiCopy += count - 1

                comp += float(present) / len(ms)
                cont += float(multiCopy) / len(ms)

            percComp = 100 * comp / len(self.markerSet)
            percCont = 100 * cont / len(self.markerSet)

        return percComp, percCont


def read_checkm_taxon_set_file(filepath):
    markers = {}
    with open(filepath) as fh:
        fh.readline()  # skip header
        # NOTE: each of the generated files only contain 2 lines...
        line = fh.readline()
    llist = line.split("\t")
    taxon_id = llist[0]
    marker_sets = read_checkm_taxon_set(line)
    for marker_set in marker_sets:
        print(
            f"{marker_set.lineageStr} lineage contains {marker_set.size()} (marker genes, marker gene sets)"
        )
        markers.update({marker_set.lineageStr: marker_set.getMarkerGenes()})
    return markers


def write_markers_keyfile(markers, fpath):
    lines = ""
    for marker in markers:
        lines += f"{marker}\n"
    with open(fpath, "w") as fh:
        fh.write(lines)


def fetch_markers_hmm(hmmfile, keyfile, outfile):
    # Now retrieve specific markers from checkm.hmm
    # For more details you may investigate CheckM source here:
    # https://github.com/Ecogenomics/CheckM/blob/e79d420194f385c67d3b3710c6beadcdf710598a/checkm/markerSets.py#L326-L343
    cmd = f"hmmfetch -f {hmmfile} {keyfile}"
    print(f"Running command: {cmd}")
    with open(outfile, "w") as stdout:
        subprocess.run(cmd, stdout=stdout, check=True, shell=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--taxon-set",
        help="File generated from <checkm taxon_set rank taxon outfile>",
        required=True
    )
    parser.add_argument(
        "--outdir",
        help="Output directory path to write taxon set lineages markers",
        required=True,
    )
    parser.add_argument(
        "--checkm-hmm",
        help="Path to checkm.hmm (Can be retrieved from https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz)",
        default="/home/evan/miniconda3/envs/sponges/checkm_data/hmms/checkm.hmm",
    )
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    markers = read_checkm_taxon_set_file(args.taxon_set)
    taxon_set_basename = os.path.splitext(os.path.basename(args.taxon_set))[0]
    
    for lineage,marker_set in markers.items():
        # First write markers keyfile
        markers_key_filepath = os.path.join(args.outdir, f"{taxon_set_basename}_{lineage}.keyfile")
        write_markers_keyfile(marker_set, markers_key_filepath)
        
        # Fetch markers hmm models using written keyfile
        markers_hmm_filepath = markers_key_filepath.replace(".keyfile", ".hmm")
        fetch_markers_hmm(
            hmmfile=args.checkm_hmm,
            keyfile=markers_key_filepath,
            outfile=markers_hmm_filepath,
        )
        print(f"Retrieved {len(marker_set)} markers hmm models in {markers_hmm_filepath}")


    all_markers = {marker for lineage, marker_set in markers.items() for marker in marker_set}
    print(f"{len(all_markers)} markers across lineages")

    # Write keyfile for hmmfetch -f <hmm> <keyfile>
    markers_key_filepath = os.path.join(args.outdir, f"{taxon_set_basename}_all_markers.keyfile")
    write_markers_keyfile(all_markers, markers_key_filepath)

    # Fetch markers hmm models using written keyfile
    markers_hmm_filepath = markers_key_filepath.replace(".keyfile", ".hmm")
    fetch_markers_hmm(
        hmmfile=args.checkm_hmm,
        keyfile=markers_key_filepath,
        outfile=markers_hmm_filepath,
    )
    print(f"Retrieved {len(all_markers)} markers hmm models in {markers_hmm_filepath}")

if __name__ == "__main__":
    main()
