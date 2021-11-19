configfile: "config/TESTconfig.yaml"

def get_inputs(wildcards):
    inputs = []
    pattern = "Results/BUSCO/CDS/{species}/{sample}/short_summary.{sample}.txt"
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
        Cor1="Results/Cor/{species}/{sample}_1.cor.fq.gz",
        Cor2="Results/Cor/{species}/{sample}_2.cor.fq.gz"
    params:
        OD="Results/Cor/{species}/"
    conda:
        "envs/rcorrector.yaml"
    shell:
        "run_rcorrector.pl -1 {input.Run1} -2 {input.Run2} -od {params.OD}"

rule removelist:
    input:
        Cor1="Results/Cor/{species}/{sample}_1.cor.fq.gz",
        Cor2="Results/Cor/{species}/{sample}_2.cor.fq.gz"
    output:
        removelist="Results/Filter/{species}/{sample}_removelist.txt"
    shell:
        """
        zgrep unfix {input.Cor1} | cut -d " " -f 1 > unfixtmp.txt
        zgrep unfix {input.Cor2} | cut -d " " -f 1 >> unfixtmp.txt
        sort -u unfixtmp.txt > {output.removelist}
        rm *tmp*
        """

rule filter:
    input:
        file="Results/Cor/{species}/{sample}_{run}.cor.fq.gz",
        removelist="Results/Filter/{species}/{sample}_removelist.txt"
    output:
        file= "Results/Filter/{species}/{sample}_{run}.sub.fq.gz"
    shell:
        "zcat {input.file} | paste -d '|' - - - - | grep -vf {input.removelist} | tr '|' '\n' | sed 's/cor$//g' | gzip > {output.file}"

rule rRNA:
    input:
        Sub1="Results/Filter/{species}/{sample}_1.sub.fq.gz",
        Sub2="Results/Filter/{species}/{sample}_2.sub.fq.gz"
    output:
        Clean1="Results/Clean/{species}/{sample}_clean.1.gz",
        Clean2="Results/Clean/{species}/{sample}_clean.2.gz",
        SAM="Results/Clean/{species}/{sample}_clean.sam",
        BAM="Results/Clean/{species}/{sample}_clean.bam"
    params:
        index="SILVA/SILVA_rRNA/SILVA_rRNA",
        met="Results/Clean/{species}/{sample}_bowtie2_metrics.txt",
        path="Results/Clean/{species}/{sample}_clean.gz"
    conda:
        "envs/mapping.yaml"
    shell:
        """
        bowtie2 --quiet --very-sensitive-local --phred33 -x {params.index} -1 {input.Sub1} -2 {input.Sub2} --threads 12 --met-file {params.met} --un-conc-gz {params.path} -S {output.SAM}
        samtools view -b {output.SAM} > {output.BAM}
        rm {output.SAM}
        """
rule trinity:
    input:
        left="Results/Clean/{species}/{sample}_clean.1.gz",
        right="Results/Clean/{species}/{sample}_clean.2.gz"
    output:
        A="Results/Trinity/{species}/{sample}/trinity_out/Trinity.fasta"
    params:
        outdir="Results/Trinity/{species}/{sample}/trinity_out"
    conda:
        "envs/trinity.yaml"
    shell:
        """
        Trinity --seqType fq --SS_lib_type RF --left {input.left} --right {input.right} --output {params.outdir} --max_memory 10G
        """
        
rule longest:
    input:
        Trinity="Results/Trinity/{species}/{sample}/trinity_out/Trinity.fasta"
    output:
        Longest="Results/Trinity/{species}/{sample}/trinity_out/Trinity.longest.fasta"
    conda:
        "envs/trinity.yaml"
    shell:
        "$TRINITY_HOME/util/misc/get_longest_isoform_seq_per_trinity_gene.pl {input.Trinity} > {output.Longest}"
        
rule TransDecoder:
    input:
        Trinity="Results/Trinity/{species}/{sample}/trinity_out/Trinity.fasta",
    output:
        "Results/Trinity/{species}/{sample}/TransDecoder/Trinity.fasta.transdecoder.cds"
    params:
        dir=directory("Results/Trinity/{species}/{sample}/TransDecoder/")
    conda:
        "envs/trinity.yaml"
    shell:
        """
        mkdir -p {params.dir}
        cd {params.dir}
        TransDecoder.LongOrfs -t ../trinity_out/Trinity.fasta
        TransDecoder.Predict -t ../trinity_out/Trinity.fasta
        """


rule Busco:
    input:
        Longest="Results/Trinity/{species}/{sample}/trinity_out/Trinity.longest.fasta",
        CDS="Results/Trinity/{species}/{sample}/TransDecoder/Trinity.fasta.transdecoder.cds"
    output:
        long="Results/BUSCO/Longest/{species}/{sample}/short_summary.{sample}.txt",
        CDS="Results/BUSCO/CDS/{species}/{sample}/short_summary.{sample}.txt"
    params:
        longout="Results/BUSCO/Longest/{sample}",
        cdsout="Results/BUSCO/CDS/{sample}",
        db="BUSCO/mammalia_odb10"
    conda:"envs/busco.yaml"
    shell:
        """
        busco -i {input.Longest} -l {params.db} -o {params.longout} -m transcriptome -c 6
        busco -i {input.CDS} -l {params.db} -o {params.cdsout} -m transcriptome -c 6
        """
        
