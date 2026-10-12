# Vilīna (विलीन) 

Finds endogenous viral elements (EVEs) in a genome assembly, annotates them, scores how
complete each one is, judges whether it could still be functional, and draws the summary
figures. Runs as a Snakemake workflow or as a single Python command.

*Vilīna*, Sanskrit for "dissolved" or "merged within": the genome's memory of ancient
infections, long since absorbed into its sequence.


## Inputs

| Input | Required | What it adds |
|---|---|---|
| `genome` (FASTA, plain or .gz) | **yes** | the search space |
| viral proteins (`viral_faa` + `viral_gpff`) | **yes** — auto-downloads NCBI RefSeq viral if left empty | hits, taxonomy (family/order/realm), genome type |
| `gtf` | no | genomic context: CDS / exon / intron / intergenic, sense/antisense, nearest gene |
| `pep` (or `cds`, translated if no pep) | no | flags annotated host genes that are themselves viral-derived (domestication evidence) |
| `cellular_decoy_faa` | no, **recommended for publication** | reciprocal search; drops EVEs whose best match is a cellular protein |
| `pirna_clusters_bed` | no | flags EVEs in piRNA clusters (non-coding functional route, important in arthropods) |

The genome alone is enough to detect EVEs and judge completeness and coding competence.
GTF and pep only add host-gene context.

## Getting the databases

```bash
bash scripts/get_refseq_viral.sh db/refseq_viral          # -> viral.protein.faa.gz + viral.protein.gpff.gz
bash scripts/get_decoy.sh swissprot db/decoy              # quick decoy (~570k proteins)
bash scripts/get_decoy.sh uniref50  db/decoy              # thorough decoy (~12 GB download)
```

Then set `viral_faa`, `viral_gpff` and `cellular_decoy_faa` in `config.yaml` to those files.

## Running

```bash
conda env create -f envs/vilina.yaml && conda activate vilina
# edit config.yaml, then either
snakemake --use-conda --cores 32 --configfile config.yaml       # parallel blastx over genome windows
# or, without Snakemake
python scripts/vilina.py run --config config.yaml
```

`envs/vilina.yaml` includes `snakemake-minimal` itself, so a fresh `conda env create` gives you
both the `snakemake` driver and everything each rule needs (`--use-conda` then just reuses this
same env per rule rather than building separate ones). Python is pinned to `3.11` on purpose —
Snakemake 7's `pulp` dependency has no build past Python 3.12, so leaving it unconstrained
(`python>=3.10`) can solve to a too-new Python that `snakemake-minimal<8` can never install
alongside, in a fresh env or an existing one.

If you created the `vilina` env **before** this file pinned Python, `snakemake: command not
found` won't fix with `conda install` into that same env — if its Python already solved to
something newer than 3.12, there is no way to add Snakemake 7 to it in place. Rebuild it instead:

```bash
conda env remove -n vilina
conda env create -f envs/vilina.yaml && conda activate vilina
```

(Install via conda/mamba, not `pip install snakemake` — pip's current `pulp` release breaks
Snakemake 7's CLI; conda-forge's metadata pins a compatible one.)

**DIAMOND ≥ 2.1.13 is required.** Version 2.1.10 reports wrong subject coordinates in
frameshift mode (`-F`), which silently breaks completeness scoring.

Test on the bundled synthetic genome (12 planted elements of known state):

```bash
python test/make_test_data.py test/data
python scripts/vilina.py run --config test/config_test.yaml   # set `diamond:` path first
python test/check_truth.py test/data/truth.tsv test/results/eves.tsv
```

## What each step does

1. **prepdb** — builds the viral DIAMOND database and a metadata table (family, order,
   realm, genome type from the RefSeq lineage; number and total length of proteins per virus).
   With a cellular decoy it also builds a combined viral + cellular database for the reverse search.
2. **split** — normalises the genome, writes a `.fai`, and cuts sequences into overlapping
   windows (default 5 Mb, 20 kb overlap) balanced across `n_chunks` blastx jobs.
3. **blastx** — `diamond blastx -F 15 --range-culling --max-hsps 0`. Frameshift-aware alignment
   lets one alignment run through the stop codons and indels typical of old EVEs, and reports them.
4. **call** — filters hits (e-value, bitscore, length; excludes retrotransposon/retrovirus families
   and cellular-homolog products by regex). Chains collinear HSPs to the same viral protein across
   gaps (insertions, introns, TE hits) into one locus. Resolves overlapping loci, keeping the best
   one and recording the others as alternative assignments. Extracts the EVE sequence and finds the
   ATG–stop ORF that best overlaps it, searching 600 bp of flank.
5. **classify** — re-aligns each ORF to the viral protein. Adds reverse search, GTF context,
   host-protein evidence and piRNA overlap, then assigns the categories below.
6. **plot / report** — figures (PNG + PDF) and a self-contained `report.html`.

## Definitions

**Completeness** (per locus, relative to the best-hit viral protein)

| class | rule |
|---|---|
| complete | ≥ 90 % of the viral protein covered **and** both N- and C-termini within 15 aa |
| near-complete | ≥ 70 % covered |
| partial | < 70 % |

`at_contig_edge` flags loci that may be truncated by the assembly rather than by biology.

**Integration blocks** (`eve_clusters.tsv`): EVEs from one viral family within 20 kb are grouped.
`genome_completeness` is the fraction of the dominant virus's proteome, in amino acids, that the
block covers. This separates "one capsid fragment" from "whole viral genome integrated".

**Coding status**

| class | rule |
|---|---|
| intact ORF | ORF ≥ 80 aa; re-aligns over ≥ 90 % of the EVE's viral span; no stops or frameshifts in the alignment |
| truncated ORF | ORF ≥ 80 aa covering ≥ 50 % of the span |
| pseudogenised | otherwise |

**Functional class** (the "potential functionality" call)

| class | meaning |
|---|---|
| likely functional (annotated host gene) | intact ORF that lies in a sense-strand annotated CDS, or whose host protein hits the same viral family — candidate domesticated gene |
| potentially functional (intact ORF) | intact ORF with no gene annotation |
| possibly functional (truncated ORF) | coding capacity partly retained |
| unlikely coding-functional | pseudogenised |

These classes measure coding *competence*, not demonstrated function. "Unlikely
coding-functional" EVEs can still act as non-coding elements, for example as piRNA-cluster
antiviral templates in arthropods. Use `in_piRNA_cluster` and small-RNA data for that. Real
evidence of function needs expression data, purifying selection (dN/dS < 1 across orthologous
EVEs in related species, e.g. with HyPhy or PAML), or conservation of the ORF across species.

**Quality flags**: `pass_reverse`/`reverse_bitscore_margin` (reverse search),
`possible_exogenous_contig` (an EVE filling a short contig may be assembled virus, not an
integration), `n_alt_assignments`/`alt_families` (ambiguous taxonomy).

## Outputs (`outdir/`)

| file | content |
|---|---|
| `eves.tsv` | one row per EVE: coordinates, taxonomy, identity, coverage, termini, stops, frameshifts, ORF, completeness, coding status, functional class, context, flags, readable `evidence` string |
| `eves.gff3`, `eves.bed` | for genome browsers; GFF3 has `viral_sequence` parents with `match_part` HSPs and `ORF` children |
| `eve_clusters.tsv` | integration blocks with viral genome completeness |
| `host_viral_homologs.tsv` | host proteins with viral hits (pep/cds input) |
| `eve_nt.fa`, `eve_orfs.faa` | EVE nucleotide sequences and ORF translations, ready for phylogenetics |
| `plots/` | 01 family × completeness · 02 family × coding potential · 03 coverage vs identity by genome type · 04 positions on chromosomes/contigs · 05 positional retention along viral proteins (N→C bias) · 06 length distribution · 07 genomic context · 08 block genome completeness |
| `report.html` | summary cards, all figures, top-200 table |
| `work/` | raw blastx, HSP table, intermediates |

## Comparative module — EVEs as rare genomic changes

A second module takes a *directory* of genomes and treats orthologous EVE integrations as
rare genomic changes (RGCs) for phylogenetics, on the same logic as SINE/LINE insertion
phylogenies: an independent germline integration at a given position is effectively a
unique, unrepeatable event, so taxa sharing one share an ancestor.

```bash
# genomes_dir: one FASTA per taxon (file stem = taxon name); set taxon_order for plots
python scripts/vilina.py compare --config config.yaml --jobs 4
# or: snakemake -s Snakefile compare --use-conda --cores 32
```

It runs the single-genome pipeline on each genome (sharing one viral database), then does
the comparative work. Use `--skip-per-taxon` to reuse per-taxon results you already have.

**The orthology problem is the whole point.** Two taxa carrying a *Bornaviridae* EVE tells
you nothing — bornaviruses integrate repeatedly and independently. What counts is the same
integration *at the same site*. So the module ignores virus identity as evidence of
homology and works from the flanks instead: it takes `flank_bp_ortho` (3 kb) either side of
every EVE, maps those flanks onto every other genome, and reconstructs where the insertion
site lies in each. Each (EVE, other taxon) pair is then scored:

| state | meaning |
|---|---|
| **1** | an EVE of the same viral family sits at the orthologous site |
| **0** | the two flanks map **adjacently** — the site exists and is demonstrably empty |
| **?** | site not locatable, or the flanks are far enough apart to hide an insertion |

The 0/? distinction is what makes the matrix defensible: only a located, demonstrably empty
site is evidence of absence. An unassembled or rearranged region scores `?`, never `0`.
Shared calls additionally have to be recovered from **both** directions
(`require_reciprocal`) before two EVEs are linked, which is what stops independent
integrations of one virus family at different sites from collapsing into a false
synapomorphy.

Linked EVEs form orthogroups; each is one binary character.

### Comparative outputs (`compare_outdir/`)

| file | content |
|---|---|
| `eve_orthogroups.tsv` | one row per orthologous insertion locus: family, member EVEs, `taxa_present`/`taxa_absent`, per-taxon state columns, `informative` flag |
| `eve_matrix.nex` | NEXUS with a PAUP block set up for **Dollo** parsimony (`ctype dollo.up:all`) |
| `eve_matrix_mrbayes.nex` | `DATATYPE=RESTRICTION` with `lset coding=noabsencesites` |
| `eve_matrix.fasta`, `eve_matrix.tsv` | binary matrix for IQ-TREE (`-st BIN -m GTR2+ASC`) or your own scripts |
| `eve_splits.tsv` | insertions supporting each bipartition, with a four-gamete compatibility check |
| `ortho_observations.tsv` | every (EVE, target taxon) site call, with coordinates and flank identities |
| `pairwise_shared.tsv` | taxon × taxon counts of shared insertions |
| `plots/C01_presence_matrix` | presence / empty / undetermined across all loci |
| `plots/C02_split_support` | insertions per bipartition, flagging those conflicting with a stronger split |
| `plots/C03_shared_insertions` | pairwise shared-insertion heatmap |
| `plots/C04_dotplot_<A>_vs_<B>` | **pairwise dotplots**: each shared EVE at its position in A vs B, coloured by viral family; open markers are located-but-empty sites |
| `plots/C05_ribbon` | **ribbon plot** in your `taxon_order`, ribbons joining orthologous insertions between adjacent taxa, coloured by viral family |
| `plots/C06_cladogram` | **cladogram** from Dollo parsimony, each branch labelled with its supporting insertions and bootstrap |
| `eve_tree.nwk` | best tree, Newick, with bootstrap as internal node labels |
| `eve_tree_all_optimal.nwk`, `eve_tree_consensus.nwk` | every equally-optimal tree, and their strict consensus (written when >1) |
| `eve_tree_support.tsv` | per clade: which insertions arose there, losses below it, bootstrap |
| `compare_report.html` | all of the above in one page |

### The cladogram

The insertions are treated as cladistic characters and the tree is built under **Dollo
parsimony**, scored purely by the number of secondary losses it forces: each character's
1-state arises exactly once, at the MRCA of the taxa carrying it. Ordinary (Fitch) parsimony
is wrong here because it is free to invent the same insertion twice — the one event the
method assumes cannot happen. `?` is a wildcard resolved toward fewest losses, so an
undetermined site never manufactures a deletion.

The search is exhaustive while the topology count allows (`exhaustive_max_trees`, default
200 000 — comfortably all topologies up to ~9 taxa), otherwise stepwise addition plus SPR
hill-climbing from several random addition orders. Outgroup monophyly is enforced. Support
comes from two independent numbers on every branch: the **count of insertions** whose single
gain maps there (bold), and **bootstrap** percentage over character resampling (grey).

```bash
# standalone, on a matrix you already have
python scripts/vilina_tree.py --groups eve_orthogroups.tsv --taxa T1,T2,T3,T4,T5 \
    --outgroup T5 --bootstrap 500 --outdir .
```

**Set `outgroup:` deliberately.** Dollo is a rooted criterion, so polarity is entirely the
root's doing, and the default (last taxon in `taxon_order`) is a placeholder, not a choice.
Rooting *inside* a group that shares insertions forces those insertions to be reconstructed
as root gains followed by losses — the tree is still right but the branch-by-branch story is
not. In the bundled test, rooting on `T5` (which belongs to the true `(T4,T5)` clade) does
exactly this: the unrooted topology still matches the true tree, but the `(T4,T5)` insertion
is reported at the root with a loss.

Read the two support numbers together. A branch carrying several reciprocally-confirmed
insertions is strong; a branch carrying one is not, and will usually tie with any conflicting
loss character — the search then returns multiple optimal trees and the consensus collapses
that node, which is the correct outcome rather than a failure.

### Reading the result

Run `eve_matrix.nex` under **Dollo parsimony**, not standard parsimony: an insertion is a
gain that happens once and can afterwards only be lost by deletion. A character that looks
like homoplasy under the four-gamete test is usually a *loss*, which Dollo handles natively
— in the bundled 5-taxon test one planted character is a secondary deletion, and it shows up
in `eve_splits.tsv` as the single "conflicting" split while the compatible splits reproduce
the true tree exactly.

Treat `C02_split_support` as the headline: it says how many independent insertions support
each branch. A handful of clean, reciprocally-confirmed insertions on a branch is strong
evidence; one is an anecdote.

Test it on the bundled 5-taxon simulation, whose true tree is `(((T1,T2),T3),(T4,T5))`:

```bash
python test/make_compare_test_data.py test/compare_data
python scripts/vilina_compare.py --config test/config_compare.yaml
python test/check_compare_truth.py test/compare_data test/results_cmp
```

It plants a root insertion, clade-defining insertions, two autapomorphies, one secondary
loss, and two *independent* same-family insertions at different sites (which must not be
merged). All nine events are recovered, the convergent pair stays separate, the compatible
splits equal the true tree's internal branches, and the cladogram's unrooted topology
matches the true tree. The checker verifies all of this, including the topology.

### Comparative caveats

- **Divergence limit.** `minimap2_preset: asm20` tops out around 5–10 % flank divergence.
  Past that, flanks stop mapping and everything scores `?`. For deeper comparisons set
  `flank_mapper: blastn` (`dc-megablast`), which is slower but far more sensitive. EVE
  insertion phylogenetics works best within a family or genus, where non-coding flanks are
  still alignable; it is not a tool for deep arthropod backbone nodes.
- **Assembly quality dominates.** Fragmented assemblies put EVEs near contig edges, losing
  one flank (presence still callable, absence not). Scaffold-level or better is strongly
  preferable, and `n_taxa_missing` in the orthogroup table tells you how much resolution
  you are losing.
- **Lineage sorting and tandem arrays.** `paralogous=1` flags groups with multiple EVEs from
  one taxon (tandem duplications of an integration). Treat those characters sceptically.
- **Reference bias.** Insertions are only discovered where a taxon has one, so taxa with
  poorer assemblies contribute fewer characters. The matrix is not a random sample.
- **Unlinked fragments.** A fragmented EVE call can yield a singleton locus that no other
  taxon could be scored against (`unlinked=1`, all other states `?`). These carry no signal
  and are excluded from `n_loci_scorable`; they are kept in the table rather than silently
  dropped so the fragmentation stays visible.
- **The tree is one line of evidence.** Insertion data is nearly homoplasy-free but sparse,
  so it resolves some nodes decisively and leaves others untouched. It complements sequence
  phylogenies at contentious nodes rather than replacing them.

## Notes and tuning

- **Sensitivity:** RefSeq viral misses much of the invertebrate RNA virosphere. For arthropods,
  add metatranscriptome-derived virus proteins as extra `viral_faa` entries with a `viral_meta`
  TSV (`acc, product, organism, family, order, realm, genome_type`). Use `--more-sensitive` for
  old, divergent EVEs.
- **Large genomes (>10 Gbp):** blastx is the main cost. Raise `n_chunks` to spread windows over
  a cluster (`snakemake --profile slurm`). `diamond_extra: "--block-size 6 --index-chunks 1"`
  speeds things up on nodes with plenty of memory. Hard-masking repeats first speeds up the
  search but can hide EVEs embedded in TEs.
- **ERVs:** Retroviridae and LTR-retrotransposon families are excluded by default because they
  swamp the output. Set `exclude_families: []` to keep them.
- **False positives:** cellular-origin genes in giant DNA viruses (kinases, ankyrins, RNR) produce
  spurious "EVEs". The product regex removes the obvious ones. The reverse search against a
  non-viral decoy (e.g. UniRef50 without Viruses) is the proper filter. Confirm key loci with
  phylogenies and by checking that host sequence flanks the integration.
