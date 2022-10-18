#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
import os, os.path

pr = argparse.ArgumentParser("Prepare workspace")
pr.add_argument("sample_annotation")
pr.add_argument("--out_dir", default=".",
    help = "output directory")
pr.add_argument("--groups", nargs="*",
    help = "grouping variables")

argv = pr.parse_args()

root_dir = os.path.abspath(argv.out_dir)


sample_annot = pd.read_csv(argv.sample_annotation, sep="\t")
sample_annot.index = sample_annot["sample_id"]


def replace_fext(fname, ext):
    return fname[:-len(ext)] + ext


# create pair annotation file

pairs = []
for name, group in sample_annot.groupby("clinical_id", sort=False):

    # find the control sample
    idx = group["sample_type"] == "Normal"
    control_ids = group.loc[idx, "sample_id"]
    if len(control_ids) < 1:
        sys.stderr.write("Missing control:\n{0}\n".format(
            list(group["sample_id"])))
        continue
    elif len(control_ids) > 1:
        sys.stderr.write(
            "Multiple controls: {0}, "
            "using first one\n".format(control_id))
    control_id = control_ids[0] 
            
    # create case-control tuples
    for i, r in group[~idx].iterrows():
        case_id = r["sample_id"]
        # if normal exists, it would have been assigned to control_id
        pairs.append((case_id, control_id))


pair_annot = pd.DataFrame(pairs, columns = ["case_id", "control_id"])
pair_annot.index = pair_annot["case_id"]
samples = list(pair_annot.index)

pair_annot.to_csv(os.path.join(root_dir, "pair_annotation.tsv"),
    sep="\t", index=False)


# create symlinks to bam files

bam_dir = os.path.join(root_dir, "bam")

if not os.path.exists(bam_dir): os.makedirs(bam_dir)
for sample in sample_annot.index:
    if pd.isnull(sample): continue
    source = sample_annot.loc[sample, "clean_bam_file_capture"]
    if pd.isnull(source): continue
    target = os.path.join(bam_dir, "{0}.bam".format(sample))
    #print("{0} -> {1}".format(source, target))
    if not os.path.exists(target):
        os.symlink(source, target)
        os.symlink(replace_fext(source, "bai"), replace_fext(target, "bai"))


# create bam lists for groups
if argv.groups:
    group_dir = os.path.join(root_dir, "group")
    if not os.path.exists(group_dir): os.makedirs(group_dir)
    for name, group in sample_annot.groupby(argv.groups):
        samples = [str(x) for x in group["sample_id"] if not pd.isnull(x)]
        paths = [os.path.join(bam_dir, "{0}.bam".format(s)) for s in samples]
        group_fname = "{0}.vtr".format('-'.join(
            [x.lower() for x in name]
        ))
        with open(os.path.join(group_dir, group_fname), 'w') as outf:
            for path in paths:
                outf.write(path + '\n')

