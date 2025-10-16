# Extract parents
awk -F'\t' 'BEGIN{OFS="\t"}
$0 !~ /^#/ && $3=="LTR_retrotransposon"{
  id=".";
  n=split($9,a,";");
  for(i=1;i<=n;i++){
    split(a[i],b,"=");
    if(b[1]=="ID"){ id=b[2]; }
  }
  print $1,$4,$5,$7,id
}' TEanno.cds.gff3 > parents.tsv

awk -F'\t' 'BEGIN{OFS="\t"}
$0 !~ /^#/ && $3=="long_terminal_repeat"{
  par="."; cid=".";
  n=split($9,a,";");
  for(i=1;i<=n;i++){
    split(a[i],b,"=");
    if(b[1]=="Parent"){ par=b[2]; }
    else if(b[1]=="ID"){ cid=b[2]; }
  }
  print $1,$4,$5,$7,par,cid
}' TEanno.cds.gff3 > ltrs.tsv

