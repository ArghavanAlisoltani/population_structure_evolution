# Extract parents
awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && $3=="LTR_retrotransposon"{
  id="."; if(match($9,/ID=([^;]+)/,m)) id=m[1];
  print $1,$4,$5,$7,id
}' TEanno.cds.gff3 > parents.tsv
# cols: seqid  start  end  strand  parent_id

# Extract child LTRs
awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && $3=="long_terminal_repeat"{
  par="."; if(match($9,/Parent=([^;]+)/,p)) par=p[1];
  id=".";  if(match($9,/ID=([^;]+)/,m)) id=m[1];
  print $1,$4,$5,$7,par,id
}' your.edta.gff3 > ltrs.tsv
# cols: seqid  start  end  strand  parent_id  ltr_id
