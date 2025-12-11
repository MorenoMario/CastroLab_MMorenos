## 1 [plasX](https://github.com/michaelkyu/PlasX)

```
anvi-gen-contigs-database -L 0 -T $THREADS --project-name $PREFIX -f unique_contigs.fasta -o $PREFIX.db
(anvio-8) mmorenos@mmorenos:012_eduCastroLab$ nohup anvi-run-ncbi-cogs   -c unique_contigs_CLAB.db  --cog-data-dir /home/mmorenos/github/anvio/anvio/data/misc/COG/COG24/COG14/ --cog-version COG14 &
```
