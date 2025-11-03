gff=split-ctg1.1Gb.fa.mod.EDTA.final/split-ctg1.1Gb.fa.mod.EDTA.intact.gff3
awk -F'\t' 'BEGIN{OFS="\t"; print "chrom","start","end","strand","elem_id","class","ltr_identity","K_JC","age_gen"}
!/^[#]/ {
  attr=$9
  # keep only LTR elements
  if (!(attr ~ /Classification=.*LTR/i)) next

  # prefer ltr_identity=; fall back to identity= if present
  pid=""
  if (match(attr, /ltr_identity=([0-9.]+)/, m)) pid=m[1]*100
  else if (match(attr, /identity=([0-9.]+)/, m)) pid=m[1]*100
  if (pid=="") next

  id="."; if (match(attr, /ID=([^;]+)/, i)) id=i[1]
  cls="LTR"; if (match(attr, /Classification=([^;]+)/, c)) cls=c[1]

  p = 1 - (pid/100.0)
  if (p >= 0.7499) next                    # JC domain guard
  K = -(3.0/4.0) * log(1 - (4.0/3.0)*p)    # JC-corrected divergence
  mu = 2.2e-9
  age = K / (2*mu)                         # generations

  print $1,$4,$5,$7,id,cls,(pid/100.0),K,age
}' "$gff" > ctg1_LTR_ages_gen.tsv

