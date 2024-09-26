# Helianthemum marifolium Bash

We are grateful to the staff of the CICA supercomputer for their
guidance in using the High Performance Computing (HPC) facility

## Cleaning:

``` bash
#!/bin/bash                 
#$ -N limpiar2                  
#$ -wd /home/amartin/Mapeo  
#son relativos a este directorio
#$ -o limpio2.salida                         
#$ -e limpio2.err                           
#$ -q corta_multicore                       
#$ -pe smp 20               
#$ -l h_vmem=4G

rm -r /scratch/amartin
mkdir -p /scratch/amartin/database
cp /home/amartin/smr_v4.3_default_db.fasta /scratch/amartin/database/

FILES=$(find . -maxdepth 1 -type f | rev | cut -c9-16 | rev | sort | uniq -d)
for i in $FILES
do
    mkdir -p /scratch/amartin/"${i}"
    cp /home/amartin/Mapeo/${i}_1.fq.gz /home/amartin/Mapeo/${i}_2.fq.gz 
    /scratch/amartin/"${i}"

    # RCORRECTOR:
    source /home/software/anaconda3/etc/profile.d/conda.sh
    conda activate Rcorrector-1.0.7
    run_rcorrector.pl -1 /scratch/amartin/"${i}"/${i}_1.fq.gz -2 
    /scratch/amartin/"${i}"/${i}_2.fq.gz -od 
    /scratch/amartin/"${i}"
    rm -r /scratch/amartin/"${i}"/${i}_1.fq.gz /scratch/amartin/"${i}"/${i}_2.fq.gz 

    # FASTP:
    conda activate fastp
    fastp -i /scratch/amartin/"${i}"/${i} -o /scratch/amartin/"${i}"/${i}_1.fq.gz -I 
    /scratch/amartin/"${i}"/${i}_2.cor.fq.gz -O /scratch/amartin/"${i}"/${i}_2.fq.gz
    --detect_adapter_for_pe --n_base_limit 5 --qualified_quality_phred 20
    --unqualified_percent_limit 50

    rm /scratch/amartin/"${i}"/${i}_1.cor.fq.gz /scratch/amartin/"${i}"/${i}_2.cor.fq.gz

    mv /home/amartin/Mapeo/fastp.html /home/amartin/Mapeo/fastp/${i}.html
    rm /home/amartin/Mapeo/fastp.json

    # SORTMERNA:
    conda activate sortmerna-4.3.6
    sortmerna --ref /scratch/amartin/database/smr_v4.3_default_db.fasta --reads ${i}_1.fq.gz
    --reads ${i}_2.fq.gz --aligned rRNA --fastx --other cleanedrRNA.fastq --threads 20 -v --out2
    --workdir /scratch/amartin/"${i}"

    mv /home/amartin/Mapeo/rRNA.log /home/amartin/Mapeo/log/${i}.log
    mv /home/amartin/Mapeo/cleanedrRNA.fastq_fwd.fq.gz /home/amartin/Mapeo/limpio/${i}_fwd.fq.gz
    mv /home/amartin/Mapeo/cleanedrRNA.fastq_rev.fq.gz /home/amartin/Mapeo/limpio/${i}_rev.fq.gz
    rm /home/amartin/Mapeo/rRNA_fwd.fq.gz  /home/amartin/Mapeo/rRNA_rev.fq.gz
    rm -r /scratch/martin/${i}
done
```

## Normalization and de novo transcriptome assembly:

### Normalization :

``` bash
insilico_read_normalization.pl  --seqType Fq  --left 
combinado_R1.fq.gz --right combinado_R2_fq.gz 
--max_cov 30 --pairs_together --PARALLEL_STATS 
--SS_lib_type FR --CPU 20
```

### Trinity:

``` bash
#$ -S /bin/bash                 
#$ -N Trinity                   
#$ -wd /home/abelardo/trinity2   
#$ -o Trinity.salida            
#$ -e Trinity.err               
#$ -q media_multicore          
#$ -pe smp 40                   
#$ -l h_vmem=4G

module load python3-3.8.2
module load trinity-2.8.6

/home/software/trinityrnaseq-2.8.6/Trinity --seqType fq --no_normalize_reads --left /home/abelardo/normalizados/noclasificado_1_2.fq.gz.normalized_K25_maxC30_minC0_maxCV10000.fq --right /home/abelardo/normalizados/noclasificado_2_2.fq.gz.normalized_K25_maxC30_minC0_maxCV10000.fq --SS_lib_type FR --max_memory 100G --output ./Trinity_output  --CPU 40
```

### Oases:

``` bash
#$ -S /bin/bash                             
#$ -N Oases                                     
#$ -wd /home/abelardo/oases            
#$ -o Oases.salida                              
#$ -e Oases.err                                 
#$ -q corta_multicore                     
#$ -pe smp 20                           
#$ -l h_vmem=4G


module load velvet-1.2.10
velveth ensamblado_oases_25  25 -fastq.gz -shortPaired -separate /home/abelardo/normalizados/noclasificado_1_2.fq.gz.normalized_K25_maxC30_minC0_maxCV10000.fq /home/abelardo/normalizados/noclasificado_2_2.fq.gz.normalized_K25_maxC30_minC0_maxCV10000.fq
velvetg ensamblado_oases_25 -read_trkg yes
/home/software/oases-0.2.8/oases ensamblado_oases_25

velveth ensamblado_oases_31  31 -fastq.gz -shortPaired -separate /home/abelardo/normalizados/noclasificado_1_2.fq.gz.normalized_K25_maxC30_minC0_maxCV10000.fq /home/abelardo/normalizados/noclasificado_2_2.fq.gz.normalized_K25_maxC30_minC0_maxCV10000.fq
velvetg ensamblado_oases_31 -read_trkg yes
/home/software/oases-0.2.8/oases ensamblado_oases_31

velveth ensamblado_oases_41  41 -fastq.gz -shortPaired -separate /home/abelardo/normalizados/noclasificado_1_2.fq.gz.normalized_K25_maxC30_minC0_maxCV10000.fq /home/abelardo/normalizados/noclasificado_2_2.fq.gz.normalized_K25_maxC30_minC0_maxCV10000.fq
velvetg ensamblado_oases_41 -read_trkg yes
/home/software/oases-0.2.8/oases ensamblado_oases_41

velveth ensamblado_oases_51  51 -fastq.gz -shortPaired -separate /home/abelardo/normalizados/noclasificado_1_2.fq.gz.normalized_K25_maxC30_minC0_maxCV10000.fq /home/abelardo/normalizados/noclasificado_2_2.fq.gz.normalized_K25_maxC30_minC0_maxCV10000.fq
velvetg ensamblado_oases_51 -read_trkg yes
/home/software/oases-0.2.8/oases ensamblado_oases_51
```

``` bash
cat trinity.fa transcripts25.fa transcripts31.fa transcripts41.fa transcripts51.fa > combinado.fa
```

## Evidentialgene:

``` bash
tr2aacds.pl -tidy -MAXMEM 60000 -log -cdna combinado.fa
```

## Quality control:

### BUSCO:

``` bash
busco -i ./combinado.okey.fa -o busco --lineage brassicales --mode transcriptome
```

### Quast:

``` bash
quast -o ./ combinado.okay.fa
```

## Mapping:

``` bash
salmon index --index index --transcripts combinado.okay.fa
```

``` bash
#!/bin/bash
FILES=$(find . -type f | rev | cut -c11-18 | rev | sort | uniq -d)
for i in $FILES
do
    samp=$(basename "$i")
    echo "Processing sample ${samp}"
  salmon quant --index index -l A -1 "${i}_fwd.fq.gz" -2 "${i}_rev.fq.gz" -p 4
  --dumpEq --hardFilter --skipQuant -o "./quants/${samp}_quant"
done
```

## Gene_level:

``` bash
gunzip -d *_quant/aux_info/eq_classes.txt.gz

corset -g 1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4 -n A1,A2,A3,A4,B1,B2,B3,B4,C1,C2,C3,C4,D1,D2,D3
-i salmon_eq_classes */aux_info/eq_classes.txt
```

## Annotation:

### Secuences coding:

``` bash
TransDecoder.LongOrfs -t combinado.okay.fa

TransDecoder.Predict -t combinado.okay.fa
```

### Homoligy and conserve domains:

``` bash
blastx -db uniprot_sprot.pep -query combinado.okay.fa -num_threads 20 -max_target_seqs 1 -outfmt 6 -evalue 1e-5 > blastx.outfmt6
```

### Transmap:

This transmap is used for annotation with Trinotate. The transcripts are
grouped according to the Corset’s clusters for annotation at the gene
level, and the transcripts eliminated by Corset due to low count numbers
are grouped in a cluster called ClusterX. This allows for their later
elimination without affecting the analysis of enrichment in GO terms.

``` bash
grep '^>' combinado.okay.fa | awk '{print $1}' > tra
sed 's/>//g' tra > tra2

cut -f1 -d$'\t' clusters.txt > tra.clu

comm -3 <(sort tra2) <(sort tra.clu) > rep

awk -F'\t' '{print $0"\tclusterx"}' rep > rep_nuevo

cat clusters.txt rep_nuevo > transmap

awk -F'\t' '{print $2 "\t" $1}' transmap  > transmapfinal
```

### Integrating:

``` bash
Trinotate --db myTrinotate.slite --init --gene_trans_map transmapfinal --transcript_fasta combinado.okay.fa --transdecoder_pep combinado.okay.fa.transdecoder.pep

Trinotate --db myTrinotate.slite --LOAD_swissprot_blastx  blastx.outfmt6

Trinotate --db myTrinotate.slite --report > trinotate_annotation.xls
```

### Getting GO Terms:

``` bash

extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls trinotate_annotation.xls --include_ancestral_terms -T > trinotate_go_annotation.txt

sed '/^clusterx/d' trinotate_annotation.xls  > trinotate_annotation.xls

extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls trinotate_annotation.xls --include_ancestral_terms -G > trinotate_go_annotation_gene.txt
```
