# Vilina (विलीन) -- endogenous viral element (EVE) discovery & characterisation workflow
#   snakemake --use-conda --cores 32 --configfile config.yaml
#
# Required: genome FASTA + viral protein DB (or set download_refseq_viral: true)
# Optional: gtf, pep (or cds), cellular_decoy_faa, pirna_clusters_bed

import os

configfile: "config.yaml"

OUT = config["outdir"]
W = str(OUT) + '/work'
DB = str(W) + '/db'
SCRIPT = os.path.join(workflow.basedir, "scripts/vilina.py")
CFG = os.path.abspath(workflow.configfiles[-1]) if workflow.configfiles else "config.yaml"
N = int(config.get("n_chunks", 8))
CHUNKS = [str(i).zfill(3) for i in range(N)]
DMD = config.get("diamond", "diamond")
SENS = config.get("diamond_sensitivity", "--sensitive")
EXTRA = config.get("diamond_extra", "") or ""
ENV = "envs/vilina.yaml"
BLASTX_FIELDS = ("qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue "
                 "bitscore slen qlen qframe qseq_translated cigar")
BLASTP_FIELDS = "qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen"


def aslist(x):
    return [] if not x else (x if isinstance(x, list) else [x])


USE_REVERSE = bool(config.get("cellular_decoy_faa"))
USE_HOST = bool(config.get("pep") or config.get("cds"))
VIRAL_FAA = aslist(config.get("viral_faa")) or [str(W) + '/refseq/viral.protein.faa.gz']
VIRAL_GPFF = aslist(config.get("viral_gpff")) or (
    [str(W) + "/refseq/viral.protein.gpff.gz"] if not config.get("viral_faa") else [])


rule all:
    input:
        str(OUT) + '/eves.tsv', str(OUT) + '/eves.gff3', str(OUT) + '/eves.bed', str(OUT) + '/eve_clusters.tsv',
        str(OUT) + '/plots/.done', str(OUT) + '/report.html'


rule download_refseq_viral:
    """Used only when viral_faa is not given in the config."""
    output:
        faa=str(W) + '/refseq/viral.protein.faa.gz',
        gpff=str(W) + '/refseq/viral.protein.gpff.gz',
    params:
        base="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral",
    shell:
        """
        cd $(dirname {output.faa})
        for f in $(curl -s {params.base}/ | grep -o 'viral\\.[0-9]*\\.protein\\.\\(faa\\|gpff\\)\\.gz' | sort -u); do
            curl -sSfO {params.base}/$f
        done
        cat viral.*.protein.faa.gz  > viral.protein.faa.gz.tmp && mv viral.protein.faa.gz.tmp viral.protein.faa.gz
        cat viral.*.protein.gpff.gz > viral.protein.gpff.gz.tmp && mv viral.protein.gpff.gz.tmp viral.protein.gpff.gz
        """


rule prepdb:
    input:
        faa=VIRAL_FAA, gpff=VIRAL_GPFF,
    output:
        dmnd=str(DB) + '/viral.dmnd', meta=str(DB) + '/viral_meta.tsv',
        rev=str(DB) + '/reverse.dmnd' if USE_REVERSE else [],
    params:
        meta='--meta ' + str(config['viral_meta']) if config.get("viral_meta") else "",
        gpff=lambda w, input: '--gpff ' + str(' '.join(input.gpff)) if input.gpff else "",
        decoy=("--decoy " + " ".join(aslist(config.get("cellular_decoy_faa")))) if USE_REVERSE else "",
    threads: workflow.cores
    conda: ENV
    shell:
        "python {SCRIPT} prepdb --faa {input.faa} {params.gpff} {params.meta} {params.decoy} "
        "--outdir {DB} --diamond {DMD} --threads {threads}"


rule split_genome:
    input: config["genome"]
    output:
        fa=str(W) + '/genome/genome.fa', fai=str(W) + '/genome/genome.fa.fai',
        chunks=expand(str(W) + '/genome/chunk_{c}.fa', c=CHUNKS),
    conda: ENV
    shell:
        "python {SCRIPT} split --genome {input} --outdir {W}/genome --n-chunks {N} "
        "--window {config[window_size]} --overlap {config[window_overlap]}"


rule blastx:
    input:
        q=str(W) + '/genome/chunk_{c}.fa', db=str(DB) + '/viral.dmnd',
    output: str(W) + '/blastx/chunk_{c}.tsv'
    threads: max(1, workflow.cores // 2)
    conda: ENV
    shell:
        """
        if [ -s {input.q} ]; then
          {DMD} blastx -q {input.q} -d {input.db} -o {output} {SENS} -F 15 --range-culling --top 10 \
            --max-hsps 0 -e {config[max_evalue]} --threads {threads} --quiet {EXTRA} \
            --outfmt 6 {BLASTX_FIELDS}
        else
          touch {output}
        fi
        """


rule call:
    input:
        hits=expand(str(W) + '/blastx/chunk_{c}.tsv', c=CHUNKS),
        meta=str(DB) + '/viral_meta.tsv', genome=str(W) + '/genome/genome.fa',
    output:
        raw=str(W) + '/eves_raw.tsv', hsps=str(W) + '/eve_hsps.tsv', nt=str(W) + '/eve_nt.fa',
        orfs=str(W) + '/eve_orfs.faa', aln=str(W) + '/eve_aligned_aa.faa',
    conda: ENV
    shell:
        "python {SCRIPT} call --hits {input.hits} --meta {input.meta} --genome {input.genome} "
        "--config {CFG} --outdir {W}"


rule orf_blastp:
    input: q=str(W) + '/eve_orfs.faa', db=str(DB) + '/viral.dmnd'
    output: str(W) + '/orf_vs_viral.tsv'
    threads: workflow.cores
    conda: ENV
    shell:
        "if [ -s {input.q} ]; then {DMD} blastp -q {input.q} -d {input.db} -o {output} {SENS} -k 25 "
        "-e 1e-5 --threads {threads} --quiet --outfmt 6 {BLASTP_FIELDS}; else touch {output}; fi"


rule reverse_blastp:
    input: q=str(W) + '/eve_aligned_aa.faa', db=str(DB) + '/reverse.dmnd'
    output: str(W) + '/reverse.tsv'
    threads: workflow.cores
    conda: ENV
    shell:
        "if [ -s {input.q} ]; then {DMD} blastp -q {input.q} -d {input.db} -o {output} {SENS} -k 25 "
        "-e 1e-3 --threads {threads} --quiet --outfmt 6 qseqid sseqid pident length evalue bitscore; "
        "else touch {output}; fi"


rule hostpep:
    input: [config[k] for k in ("pep", "cds") if config.get(k)][:1]
    output: str(W) + '/host_pep.faa'
    params:
        arg=('--pep ' + str(config['pep']) if config.get("pep") else '--cds ' + str(config.get('cds'))),
    conda: ENV
    shell: "python {SCRIPT} hostpep {params.arg} --out {output}"


rule host_blastp:
    input: q=str(W) + '/host_pep.faa', db=str(DB) + '/viral.dmnd'
    output: str(W) + '/host_vs_viral.tsv'
    threads: workflow.cores
    conda: ENV
    shell:
        "if [ -s {input.q} ]; then {DMD} blastp -q {input.q} -d {input.db} -o {output} {SENS} -k 25 "
        "-e 1e-5 --threads {threads} --quiet --outfmt 6 {BLASTP_FIELDS}; else touch {output}; fi"


rule classify:
    input:
        raw=str(W) + '/eves_raw.tsv', meta=str(DB) + '/viral_meta.tsv', orf=str(W) + '/orf_vs_viral.tsv',
        rev=str(W) + '/reverse.tsv' if USE_REVERSE else [],
        host=str(W) + '/host_vs_viral.tsv' if USE_HOST else [],
    output:
        tsv=str(OUT) + '/eves.tsv', gff=str(OUT) + '/eves.gff3', bed=str(OUT) + '/eves.bed',
        clus=str(OUT) + '/eve_clusters.tsv', hv=str(OUT) + '/host_viral_homologs.tsv', summ=str(OUT) + '/summary.json',
    params:
        extra=" ".join(filter(None, [
            '--reverse-hits ' + str(W) + '/reverse.tsv' if USE_REVERSE else "",
            '--host-hits ' + str(W) + '/host_vs_viral.tsv' if USE_HOST else "",
            ("--gtf " + str(config["gtf"])) if config.get("gtf") else "",
            '--pirna ' + str(config['pirna_clusters_bed']) if config.get("pirna_clusters_bed") else "",
        ])),
    conda: ENV
    shell:
        "python {SCRIPT} classify --indir {W} --meta {input.meta} --orf-hits {input.orf} "
        "--config {CFG} --outdir {OUT} {params.extra} && cp {W}/eve_nt.fa {W}/eve_orfs.faa {OUT}/"


rule plot:
    input: eves=str(OUT) + '/eves.tsv', clus=str(OUT) + '/eve_clusters.tsv'
    output: touch(str(OUT) + '/plots/.done')
    conda: ENV
    shell: "python {SCRIPT} plot --eves {input.eves} --clusters {input.clus} --outdir {OUT}/plots"


rule report:
    input: eves=str(OUT) + '/eves.tsv', summ=str(OUT) + '/summary.json', plots=str(OUT) + '/plots/.done'
    output: str(OUT) + '/report.html'
    params: title=config.get("sample", "genome")
    conda: ENV
    shell:
        "python {SCRIPT} report --eves {input.eves} --summary {input.summ} --plots {OUT}/plots "
        "--title '{params.title}' --out {output}"


# ---------------------------------------------------------------------------
# Comparative module (EVEs as rare genomic changes). Needs `genomes_dir` in the
# config. Per-taxon EVE calling is driven inside the module so it can share one
# viral database; run it with:   snakemake -s Snakefile compare --cores 32
# ---------------------------------------------------------------------------
COMPARE_OUT = config.get("compare_outdir") or (str(OUT).rstrip("/") + "_compare")
COMPARE_SCRIPT = os.path.join(workflow.basedir, "scripts/vilina_compare.py")


rule compare:
    """Multi-genome EVE orthology, character matrix and comparative plots."""
    output:
        groups=str(COMPARE_OUT) + '/eve_orthogroups.tsv',
        nex=str(COMPARE_OUT) + '/eve_matrix.nex',
        splits=str(COMPARE_OUT) + '/eve_splits.tsv',
        report=str(COMPARE_OUT) + '/compare_report.html',
    threads: workflow.cores
    conda: ENV
    shell:
        "python {COMPARE_SCRIPT} --config {CFG} --jobs 1"
