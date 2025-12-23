python plot_tmrca_bin_bars_v3.py \
  --input tmrca_fast_bins_full_v2_mybins_v6/tmrca_bins_frequencies.tsv \
  --out-prefix tmrca_fast_bins_full_v2_mybins_v6/tmrca_bins_barplots_v3_py_ \
  \
  --col-bin bin_label \
  --col-snp-frac frac_snps \
  --col-seg-frac frac_segments \
  --col-gene-frac frac_genes \
  --col-te-frac frac_TEs_all \
  --col-copia-frac frac_TEs_Copia \
  --col-gypsy-frac frac_TEs_Gypsy \
  \
  --x-gap 0.02 \
  --col-snp "grey" --col-gene "blue" \
  --col-te "red" --col-copia "green" --col-gypsy "red" \
  --alpha-snp 0.5 --alpha-gene 0.5 --alpha-te 0.5 \
  --alpha-copia 0.5 --alpha-gypsy 0.5 \
  --width-snp 0.2 --width-seg 0.2 --width-gene 0.4 \
  --width-te 0.6 --width-copia 0.6 --width-gypsy 0.8 \
  \
  --ybreak-main "0.05,0.20" \
  --make-area \
  --smooth-area-window 3 \
  --legend-loc right-out
