You can build a [halo](https://github.com/gear-genomics/halo) singularity container (SIF file) using

`sudo singularity build halo.sif halo.def`

Once you have built the container you can run analysis using

`singularity exec halo.sif halo count -g ref.fa sc1.bam sc2.bam`
