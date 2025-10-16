# GFF3 columns: seqid source feature start end score strand phase attributes
# Pull seqid, start, end, feature, strand, and TE Name= (family ID)
awk 'BEGIN{OFS="\t"} !/^#/{
  name="NA";
  if (match($9, /(^|;)Name=([^;]+)/, a)) name=a[2];
  print $1,$4,$5,$3,$7,name
}' your.edta.gff3 > gff_instances_family.tsv
