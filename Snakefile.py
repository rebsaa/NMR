configfile: "config/config.yaml"

def get_inputs(wildcards):
    inputs = []
    pattern = "Results/Trinity/{species}/{sample}/Diamond/out_Heterocephalus_glaber_male.tsv"
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
        PEP="Results/Trinity/{species}/{sample}/TransDecoder/Trinity.fasta.transdecoder.pep",
        CDS="Results/Trinity/{species}/{sample}/TransDecoder/Trinity.fasta.transdecoder.cds"
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
        long="Results/BUSCO/Longest/{species}/{sample}/short_summary.specific.mammalia_odb10.{sample}.txt",
        CDS="Results/BUSCO/CDS/{species}/{sample}/short_summary.specific.mammalia_odb10.{sample}.txt"
    params:
        longout="Results/BUSCO/Longest/{species}",
        cdsout="Results/BUSCO/CDS/{species}",
        db="BUSCO/mammalia_odb10"
    conda:"envs/busco.yaml"
    shell:
        """
        busco -i {input.Longest} -l {params.db} --out_path {params.longout} -o {wildcards.sample} -m transcriptome -c 6 -f
        busco -i {input.CDS} -l {params.db} --out_path {params.cdsout} -o {wildcards.sample} -m transcriptome -c 6 -f
        """

rule Diamond: 
    input:
        PEP="Results/Trinity/{species}/{sample}/TransDecoder/Trinity.fasta.transdecoder.cds"
    output:
        "Results/Trinity/{species}/{sample}/Diamond/out_Heterocephalus_glaber_male.tsv"
    params:
        INFILE="Results/Trinity/{species}/{sample}/Diamond/Diamond_input.cds",
        TMP="Results/Trinity/{species}/{sample}/Diamond/Diamond_input.tmp",
        DIR="Results/Trinity/{species}/{sample}/Diamond/"
    conda:"envs/diamond.yaml"
    shell:
        """
        cp {input.PEP} {params.INFILE}
        for DB in $(cut -d '/' -f 9 Diamond/diamond_ref_pep.txt| cut -d '.' -f 1); do
        diamond blastx -q {params.INFILE} -o {params.DIR}out_$DB.tsv -d Diamond/$DB.dmnd -p 12 --max-target-seqs 1 --id 90
        cp {params.INFILE} {params.TMP}
        cat {params.TMP} | tr '\\n' '|' | sed 's/>/\\n>/g'| grep -wvf <(cut -f 1 {params.DIR}*.tsv) | tr '|' '\\n' | sed -r '/^\\s*$/d' > {params.INFILE}
        done
        
        for DB in $(cut -d '/' -f 9 Diamond/diamond_ref_pep.txt| cut -d '.' -f 1); do
        diamond blastx -q {params.INFILE} -o {params.DIR}out_$DB.tsv -d Diamond/$DB.dmnd -p 12 --max-target-seqs 1 --id 80
        cp {params.INFILE} {params.TMP}
        cat {params.TMP} | tr '\\n' '|' | sed 's/>/\\n>/g'| grep -wvf <(cut -f 1 {params.DIR}*.tsv) | tr '|' '\\n' | sed -r '/^\\s*$/d' > {params.INFILE}
        done
        """
