```
awk '$1~/ENSMUSG/' lung_expression.txt |
    perl test_mus.pl \
        >Mmu_lung_expr.tsv

pigz -dc gencode.vM25.pc_translations.fa.gz |
    perl test_mus1.pl \
        Mmu_lung_prom.tsv |
    awk '$3>=0' >Mmu_lung_ENSP.tsv

perl test_mus2.pl \
    Mmu_lung_expr.tsv \
    <Mmu_lung_ENSP.tsv \
    >Mmu_lung_expr_prom.tsv
```