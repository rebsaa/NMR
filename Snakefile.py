configfile: "config/config2.yaml"

def get_inputs(wildcards):
    inputs = []
    pattern = "Trinity/{species}/{sample}/trinity_out/Trinity.fasta"
    for SMP in config["samples"]:
        species, sample, run = SMP.split("/")
        inputs.append(pattern.format(species=species, SMP=SMP, sample=sample, run=run))
    return inputs

rule all:
    input:
        get_inputs

rule rcorrector:
    input:
        Run1="{species}/{sample}_1.fastq.gz",
        Run2="{species}/{sample}_2.fastq.gz"
    output:
        Cor1="Cor/{species}/{sample}_1.cor.fq.gz",
        Cor2="Cor/{species}/{sample}_2.cor.fq.gz"
    conda:"envs/rcorrector.yaml"
    shell:
        "perl ./.snakemake/conda/f39c077a1e56bf61669e242ef058cbf5/bin/run_rcorrector.pl -1 {input.Run1} -2 {input.Run2} -od /data/scratch/saager94/NMR/Cor/{wildcards.species}/"

rule removelist:
    input:
        Cor1="Cor/{species}/{sample}_1.cor.fq.gz",
        Cor2="Cor/{species}/{sample}_2.cor.fq.gz"
    output:
        removelist="Filter/{species}/{sample}_removelist.txt"
    shell:
        """
        zgrep unfix {input.Cor1} | cut -d " " -f 1 > unfixtmp.txt
        zgrep unfix {input.Cor2} | cut -d " " -f 1 >> unfixtmp.txt
        sort -u unfixtmp.txt > {output.removelist}
        rm *tmp*
        """

rule filter:
    input:
        file="Cor/{species}/{sample}_{run}.cor.fq.gz",
        removelist="Filter/{species}/{sample}_removelist.txt"
    output:
        file= "Filter/{species}/{sample}_{run}.sub.fq.gz"
    shell:
        "zcat {input.file} | paste -d '|' - - - - | grep -vf {input.removelist} | tr '|' '\n' | sed 's/cor$//g' | gzip > {output.file}"

rule rRNA:
    input:
        Sub1="Filter/{species}/{sample}_1.sub.fq.gz",
        Sub2="Filter/{species}/{sample}_2.sub.fq.gz"
    output:
        Clean1="Clean/{species}/{sample}_clean.1.gz",
        Clean2="Clean/{species}/{sample}_clean.2.gz"
    params:
        index="SILVA/SILVA_rRNA"
    conda:"envs/mapping.yaml"
    shell:
        "bowtie2 --quiet --very-sensitive-local --phred33 -x {params.index} -1 {input.Sub1} -2 {input.Sub2} --threads 12 --met-file Clean/{wildcards.species}/{wildcards.sample}_bowtie2_metrics.txt --un-conc-gz Clean/{wildcards.species}/{wildcards.sample}_clean.gz -S Clean/{wildcards.species}/{wildcards.sample}_clean.sam"

rule trinity:
    input:
        left="Clean/{species}/{sample}_clean.1.gz",
        right="Clean/{species}/{sample}_clean.2.gz"
    output:
        A="Trinity/{species}/{sample}/trinity_out/Trinity.fasta"
    conda:"envs/trinity.yaml"
    shell:
        """
        Trinity --seqType fq --SS_lib_type RF --left {input.left} --right {input.right} --output Trinity/{wildcards.species}/{wildcards.sample}/trinity_out --max_memory 10G
        """
        
rule longest:
    input:
        Trinity="Trinity/{species}/{sample}/trinity_out/Trinity.fasta"
    output:
        Longest="Trinity/{species}/{sample}/trinity_out/Trinity.longest.fasta"
    conda:"envs/trinity.yaml"
    shell:
        "
        get_longest_isoform_seq_per_trinity_gene.pl {input.Trinity} > {output.Longest}
        "

rule TransDecoder:
    input:
        Trinity="Trinity/{species}/{sample}/trinity_out/Trinity.fasta"
    output:"Trinity/{species}/{sample}/TransDecoder/Trinity.fasta.transdecoder.cds"
    conda:"envs/trinity.yaml"
    shell:
        "
        TransDecoder.LongOrfs -t {input.Trinity} --output_dir TransDecoder
        TransDecoder.Predict -t {input.Trinity}  --output_dir TransDecoder
        "
        
rule Busco:
    input:
        Longest="Trinity/{species}/{sample}/trinity_out/Trinity.longest.fasta"
        CDS="Trinity/{species}/{sample}/TransDecoder/Trinityx.fasta.transdecoder.cds"
    output:
        long="BUSCO/Longest/{sample}/short_summary.{sample}.txt"
        CDS0"BUSCO/CDS/{sample}/short_summary.{sample}.txt"
    params:
        longout="BUSCO/Longest/{sample}"
        cdsout="BUSCO/CDS/{sample}"
        db="BUSCO/mammalia_odb10"
    conda:"envs/busco.yaml"
    shell:
        "
        busco -i {input.Longest} -l {params.db} -o {params.longout} -m transcriptome -c 6
        busco -i {input.CDS} -l ./mammalia_odb10 -o {params.cdsout} -m transcriptome -c 6
        "
