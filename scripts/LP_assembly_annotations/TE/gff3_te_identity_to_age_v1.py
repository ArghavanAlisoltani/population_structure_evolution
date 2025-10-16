#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert EDTA GFF3 `identity` to a rough TE insertion-age estimate.

Given a GFF3 with features annotated by EDTA and an `identity=<float>` attribute
(identity to consensus), we compute:
  divergence D = max(0, 1 - identity)
  age (generations) = D / (2 * mu)
  age (years)       = age_gens * generation_time_years

The rate of nucleotide substitution (r) used for gymnosperm species was 2.2 × 10–9

This is a *rough proxy*. For proper ages of LTR retrotransposons, use paired LTR
5'/3' divergence when available. Here we only have top-level features with
`identity`, so we provide this approximation.

Usage
-----
python gff3_te_identity_to_age_v1.py \
  --gff ../Lpine.gff3 \
  --mu 2.2e-9 \
  --generation-years 1 \
  --include all \
  --out te_age_from_identity.tsv

Arguments
---------
--gff : EDTA GFF3 file
--mu  : Neutral substitution rate per site per generation (default 1.3e-8)
--generation-years : Years per generation (default 1)
--include : one of {all,ltr,with_identity} (default ltr)
            - all: all rows regardless of ontology
            - ltr: only rows with ontology containing 'LTR' (e.g., 'LTR_retrotransposon', 'Copia_LTR_retrotransposon')
            - with_identity: any row that has an identity attribute
--out : Output TSV path (default stdout)

Output columns
--------------
seqid, start, end, strand, feature, classification, ontology, id, name,
identity, divergence, age_generations, age_years

Notes
-----
- We cap identity to [0,1]. Rows without identity are dropped unless --include=all.
- Coordinates are passed through as in the GFF3.
"""
from __future__ import annotations
import sys, argparse, re, math, csv

IDENT_RE = re.compile(r"(?:^|;)identity=([^;]+)")
CLASS_RE = re.compile(r"(?:^|;)classification=([^;]+)")
ID_RE    = re.compile(r"(?:^|;)ID=([^;]+)")
NAME_RE  = re.compile(r"(?:^|;)Name=([^;]+)")


def parse_gff(gff_path: str):
    with open(gff_path, 'r') as f:
        for line in f:
            if not line or line.startswith('#'):
                # skip comments/header but keep if it's a pragma line
                if line.startswith('##FASTA'):
                    break
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attrs = parts
            yield {
                'seqid': seqid,
                'source': source,
                'feature': feature,
                'start': int(float(start)),
                'end': int(float(end)),
                'score': score,
                'strand': strand,
                'phase': phase,
                'attrs': attrs,
            }


def get_attr(re_pat, s, default=None):
    m = re_pat.search(s)
    return m.group(1) if m else default


def clamp01(x: float) -> float:
    if math.isnan(x):
        return x
    return 0.0 if x < 0 else 1.0 if x > 1 else x


def run(args):
    rows = []
    for rec in parse_gff(args.gff):
        ontology = rec['feature']  # EDTA puts ontology in the 3rd column in this export
        has_ltr = ('LTR' in ontology)
        identity_s = get_attr(IDENT_RE, rec['attrs'])
        has_identity = identity_s is not None
        if args.include == 'ltr' and not has_ltr:
            continue
        if args.include == 'with_identity' and not has_identity:
            continue
        identity = float(identity_s) if has_identity else float('nan')
        identity = clamp01(identity) if math.isfinite(identity) else identity
        divergence = (1.0 - identity) if math.isfinite(identity) else float('nan')
        age_gens = (divergence / (2.0 * args.mu)) if (math.isfinite(divergence) and args.mu>0) else float('nan')
        age_years = (age_gens * args.generation_years) if math.isfinite(age_gens) else float('nan')
        classification = get_attr(CLASS_RE, rec['attrs'], '')
        rid = get_attr(ID_RE, rec['attrs'], '')
        name = get_attr(NAME_RE, rec['attrs'], '')
        rows.append({
            'seqid': rec['seqid'],
            'start': rec['start'],
            'end': rec['end'],
            'strand': rec['strand'],
            'feature': ontology,
            'classification': classification,
            'ontology': ontology,
            'id': rid,
            'name': name,
            'identity': identity if math.isfinite(identity) else '',
            'divergence': divergence if math.isfinite(divergence) else '',
            'age_generations': age_gens if math.isfinite(age_gens) else '',
            'age_years': age_years if math.isfinite(age_years) else '',
        })

    out_fields = ['seqid','start','end','strand','feature','classification','ontology','id','name',
                  'identity','divergence','age_generations','age_years']
    if args.out:
        with open(args.out, 'w', newline='') as w:
            writer = csv.DictWriter(w, fieldnames=out_fields, delimiter='\t')
            writer.writeheader()
            for r in rows:
                writer.writerow(r)
        print(f"[ok] wrote: {args.out}")
    else:
        writer = csv.DictWriter(sys.stdout, fieldnames=out_fields, delimiter='\t')
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description='Convert EDTA GFF3 identity to rough insertion age')
    p.add_argument('--gff', required=True, help='EDTA GFF3 file path')
    p.add_argument('--mu', type=float, default=1.3e-8, help='Substitution rate per site per generation (default 1.3e-8)')
    p.add_argument('--generation-years', type=float, default=1.0, help='Years per generation (default 1)')
    p.add_argument('--include', choices=['all','ltr','with_identity'], default='ltr', help='Which records to include (default: ltr)')
    p.add_argument('--out', help='Output TSV (default: stdout)')
    return p.parse_args(argv)

if __name__ == '__main__':
    run(parse_args())


