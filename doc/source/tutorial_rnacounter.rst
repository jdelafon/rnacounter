Using rnacounter to count reads in genomic intervals
====================================================

`Rnacounter` estimates abundances of genes and their different transcripts
from read alignments. Exons and introns can also be quantified.

It provides fast read counting in annotated genomic features as well as a simple,
yet efficient solution to the quantification of isoforms from RNA-seq data.
The method used is described in [<ref>].
A typical run is expected to take less than 2 minutes for a 1Gb BAM file from mouse
RNA sequencing, increasing linearly with the BAM size.

For all these tasks it only requires a BAM file from a read mapping on the genome,
and a single GTF/GFF file describing the exon structure
such as those provided by Ensembl or GenRep.

It is not meant to be used as a library, but through its command-line tool "rnacounter".

The code project is hosted in Github (https://github.com/delafont/rnacounter), GPL-2 licensed.

Basic Usage
-----------
The GTF is assumed to be sorted at least w.r.t. chromosome name,
and the chromosome identifiers in the GTF must be the same as the BAM references.
The BAM is supposed to be sorted and indexed. If the index is not found it will be
created automatically (takes a few minutes).

The following command will create a tab-delimited text file of gene counts,
RPKM and annotation information such as the genomic location::

   rnacounter test.bam test.gtf > counts_table.txt

Many options are then available, listed in the Options section below.

The special command `rnacounter test` runs the program on already included
small sample files to check that it works correctly::

   rnacounter test

Use `rnacounter join` to merge several output files produced using **the same annotation**,
to create a single table with counts from all samples::

   rnacounter join tab1.txt tab2.txt ... > tab_all.txt

The tables to join must have exactly the same number of lines -
any kind of filtering can break this requirement.

You can print this help page directly with::

   rnacounter -h

Input files format
------------------
A BAM file is the result of any common read aligner (Bowtie, BWA, ...).
BAM files specification:

http://samtools.github.io/hts-specs/SAMv1.pdf

A GTF file (sometimes also called GFF or GFF2/3) looks like this::

    chr1	Ensembl	exon	11868	12227	.	+	.	exon_id "ENSE00002234944"; transcript_id "ENST00000456328"; gene_id "ENSG00000223972"; gene_name "DDX11L1"
    chr1	Ensembl	exon	12009	12057	.	+	.	exon_id "ENSE00001948541"; transcript_id "ENST00000450305"; gene_id "ENSG00000223972"; gene_name "DDX11L1"
    ...

The last column contains various optional annotation, in an arbitrary order.
I exon_id is not present, a random one will be generated.
If any of the others mentioned in the example are
not present (e.g. transcript_id), they will be given the value of exon_id.
All other annotations will be ignored.
GTF/GFF specification:

http://www.ensembl.org/info/website/upload/gff.html

One can dowload suitable GTF files from Ensembl (under "Gene sets"):

http://www.ensembl.org/info/data/ftp/index.html

Output format
-------------

The output is a tab-delimited text file with the following fields:

* ID : unique feature ID, as in the gene_id/exon_id/transcript_id GTF fields.
* Count : raw read count, without any normalization.
* RPKM : count divided by the region's size and divided by the total number of reads
  (unless :option:`-n` is specified).
* Chrom, Start, End, Strand : region localization in the genome.
* GeneName : gene symbol (may not be unique).
* Type : gene, exon, intron or transcript, in case they are mixed in the output (see `-t`).
* Sense : for strand-specific protocols (see :option:`-s`), sense and antisense read counts
  are reported on two different lines. If not specified, returns '.' .
* Synonym : list of synonym IDs (see :option:`--exon_cutoff`), if any, '.' otherwise.

The file starts with a 1-line header listing the field names.

Options
-------

* :option:`-h`, :option:`--help` and :option:`-v`, :option:`--version`:

  Display information about the program usage / the version currently installed.

* :option:`-s`, :option:`--stranded`:

  If the protocol was strand-specific and this option is provided,
  sense and antisense counts are both reported in two consecutive lines
  with a different tag in the last column.
  They can be split afterwards by piping the result as for instance with
  `... | grep 'antisense'`.
  Using the `--threshold` option together with `--stranded`
  will exclude only elements with both sense and antisense counts under the threshold.

* :option:`-n`, :option:`--normalize`:

  RPKM are automatically calculated together with raw read counts. RPKM are counts
  divided by the length of the transcript as well as by a sample-specific
  normalization constant, usually the total number of aligned reads in the sample (default).
  This value can be changed to a user-defined integer.
  Typically, if you want to compare the same gene in several samples,
  the normalization will cancel out anyway
  and giving `-n 1` will speed up the process since it will skip counting the alignments.
  Some stats programs also require raw counts anyway and do their own normalization.
  To get FPKM instead, see `--fraglength`.

* :option:`-f`, :option:`--fraglength`:

  Since in a transcript of length L there are only L-F+1 different positions where
  a fragment of length F can be cut, one may want to correct for this bias before RPKM
  calculation (then usually called FPKM). Typical fragment lengths are around 350nt;
  default value is 1 (no correction). This is not to be confused with the read length.
  This option can be applied only at the gene- or transcript level.

* :option:`--nh`:

  A flag "NH" can be added to BAM files to indicate the number of times the read
  could be mapped to different locations in the genome. Adding this option
  will take this number into account by adding 1/NH instead of 1 to an exon read count.

* :option:`--noheader`:

  By default the program adds one line with column descriptors on top of the output file.
  For easier piping the result to some other program, one can choose
  not to add the header by adding this option.

* :option:`--exon_cutoff`:

  Often the annotation contains (sometimes artificial) transcript structures that are
  very close to each other and are thus hard to dinstinguish for any model due to
  the read length constraint and lack of coverage on small regions, reducing
  the model's power.
  To address this, one can merge transcripts differing by exonic
  regions of less than that many nucleotides.
  In the output, all similar transcripts are reported with the same score,
  in succession, but synonyms are listed in a supplementary column.
  Synonyms include the feature itself, in order to easily group synonym features.
  The duplicate scores are not accounted for in any calculation.
  A zero cutoff (default) disables transcripts filtering, and is especially suitable to
  "local" alignments, or to a bigger number to reduce the transcripts variety.
  A negative value sets the cutoff to read length.

* :option:`--threshold`:

  Features with counts inferior or equal to the given threshold (positive number)
  will not be reported in the ouput. By default everything is reported
  - even with zero counts.

* :option:`--gtf_type`:

  Usually one uses standard (Ensembl etc.) GTF files to count reads in
  exons/genes/transcripts. The only lines of interest are then the ones with
  value "exon" (default) in the 3rd column. If you are counting something else
  or provided your own, differently formatted GTF, with this option you can specify
  the 3rd column value of the lines to consider.

* :option:`--format`:

  One can also give an annotation file in BED format with 4 fields
  (chromosone, start, end, name), in which case each line
  is considered as an independant, disjoint interval with no splicing structure.
  Default is "gtf", can be changed to "bed".
  The 4th column of the BED format (name) must contain *unique* IDs.
  If the input format is "bed", the program cannot know which type of intervals
  is represented, thus will always report them as 'genes' in the output.
  Consistently, it cannot be used in conjunction with the :option:`--type` option.
  Since every interval in BED format is treated independently, this mode is usually
  slower (no clever features grouping).

* :option:`-t`, :option:`--type`:

  The type of feature you want to count reads in. Can be "genes" (default),
  "transcripts", "exons", "exon_frags" or "introns".
  "exon_frags" means "disjoint exon fragments", opposed to "exons" which can overlap
  (see the Overlaps,... section below).
  One can give multiple comma-separated values, in which case all
  the different features will be mixed in the output but can easily be split
  using the last column tag, as for instance with `... | grep 'exon'`.
  Then if :option:`--method` is specified it must have the same number of values as
  :option:`--type`,
  also as a comma-separated list, or a single one that is applied to all types.

* :option:`-c`, :option:`--chromosomes`:

  Consider only a subset of all chromosomes by providing a comma-separated list
  of chromosome names (that must match those of the GTF and BAM).

* :option:`-o`, :option:`--output`:

  The output is `stdout` by default (output directly to screen unless redirected).
  Alternatively one can redirect the standard output to
  a file using this option. If the file name already exists, it will be overwritten.

* :option:`-m`, :option:`--method`:

  Feature counts are inferred from the counts on (slices of) exons
  with the chosen :option:`--method`: "raw" (htseq-count-like) or
  "nnls" (non-negative least squares, see [<ref>]).
  The default is "raw" to not disturb habits, but "nnls" is advised
  especially at the transcripts level (see the Overlaps,... section below).
  For genes (:option:`-t genes`), a special method "indirect-nnls" exists that
  calculates transcripts expressions by NNLS and returns the gene count as the
  sum of all its transcripts counts.

Overlaps, redundancies and choice of counting method
----------------------------------------------------

Viewed as annotated segments along the genome line, exons often overlap
with each other, and so do the transcripts they constitute, hence the
deconvolution problem that we solve by the NNLS method (see `--method`).
The `-t` option allows to count in "genes", "transcripts", "exons",
"exon_frags" or "introns".

In "raw" mode (total number of reads aligned to the genomic region),
"exons" counts will be redundant where overlaps occur, as their reads
are counted twice. For that reason we first decompose exons in disjoint
slices "exon_frags", so that the sum of all slices counts for a gene is the
total number of reads aligned to that gene.
For the same reason, in "raw" mode, "transcripts" counts will be redundant.

Since the overlapping regions are usually big compared to the size of the exons,
discarding all ambiguous reads would remove almost all the information.
The problem does not really apply to overlaps between gene annotations,
which are usually small. Thus in "raw" mode, intersections between genes
are removed before counting. If the data is strand-specific (see `--stranded`),
there is no more ambiguity and nothing is removed.

Now for "exons", "transcripts" and not strand-specific "genes",
using the "nnls" counting mode (see :option:`--method`) will remove the redundancy
(adding some error proportional to the ambiguity, but still less inaccurate than
redundant counts).
Try::

    rnacounter test -t transcripts

compared to::

    rnacounter test -t transcripts -m nnls

The "nnls" method is equivalent to "raw" for disjoint intervals such as
"introns" and "exon_frags".

Counting in genes is traditionally done as via the "raw" method.
However, to remain consistent the expression of the gene should be
the quantity of RNA transcribed from this gene. Thus it probably makes sense
to first calculate transcripts expressions by NNLS and sum them to obtain the
gene count, which is implemented in the "indirect-nnls" method.

In summary, one should:

* count genes with either "raw" or "nnls" or "indirect-nnls", depending on one's beliefs;
* count transcripts with "nnls";
* count exons with "nnls", or rather consider the "exon_frags" instead;
* and the method of choice does not matter for introns and "exon_frags".

Miscellaneous notes
-------------------

* Multiple alignments:

  Rather than an option/default to remove multiply mapping reads, this filtering
  - if desired - should be done at the mapping step choosing the right parameters,
  or the BAM file can be filtered afterwards. On the contrary if you want to keep
  multiple mapping but correct for it, you can use the `--nh` option.

* Exons and introns annotation:

  If no exon_id is present in the GTF, a random, unique one is assigned.
  "exon_frags" names in the output table are formatted as
  "exon1|exon2" if the fragment is spanned by exon1 and exon2.

  Intronic regions that are also annotated as exons in some alternative transcripts are
  ignored whatever the chosen method is - i.e. only regions that are intronic in
  all alternative transcripts are reported.
  Because they don't have official IDs, introns slices are given names following
  this pattern: "<n>I-<gene_id>", if it is the n-th intron of that gene.
  Report to their coordinates to identify them more precisely.

* Non-integer counts:

  The fact that some reads cross exon boundaries as well as considering the NH flag
  make the reported numbers not be integers. Some discrete distributions-based
  programs for differential expression analysis require to round them.

* Custom input:

  If your GTF does not represent exons but custom genomic intervals to simply count
  reads in, provide at least a unique `exon_id` in the attributes as a feature name,
  and the type field (column 3) must be set to 'exon' or specified with the
  `--gtf_ftype` option. If not specified, `gene_id`, `transcript_id` and `exon_id`
  will all get the value of `exon_id`.

* Paired-end support:

  At the moment alignments of paired-end reads are not treated specially, i.e.
  all reads are considered as single-end.

Examples
--------

* Probably the best way to get isoforms counts::

    rnacounter -t transcripts -m nnls --nh -f 350 sample.bam mouse.gtf > transcript_counts.txt

* Compare gene counts between two conditions, HTSeq-like::

    rnacounter group1.bam mouse.gtf > gene_counts1.txt
    rnacounter group2.bam mouse.gtf > gene_counts2.txt
    rnacounter join gene_counts1.txt gene_counts2.txt > gene_counts.txt

  Then send it to DESeq/EdgeR/whatever other stats program that asks for such a table.

FAQ & Troubleshooting
---------------------

Any bug report, usage issue or feature request not listed below can be addressed to
julien.delafontaine@epfl.ch or webmaster.bbcf@epfl.ch .

* The program ends without an error but the output file is empty:

  Most probably there is a mismatch between the BAM and the annotation files,
  usually not using the same assembly, or not referencing the same chromosome names.

* I don't get the same numbers as with htseq-count:

  First check if the input data is strand-specific (htseq-count has the -s=yes by default).
  Secondly, rnacounter does not discard reads crossing exon boundaries - but adds a fraction
  of the read to the nucleotide count.
  In NNLS mode, though, the counts are expected to be significantly different
  in regions where exons overlap, since it does not remove/ignore the overlaps.

Reference
---------

<?>

