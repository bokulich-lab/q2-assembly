# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.core.type import Bool, Choices, Float, Int, List, Range, Str

megahit_params = {
    "presets": Str,
    "min_count": Int % Range(1, None),
    "k_list": List[Int % Range(15, 255, inclusive_end=True)],
    "k_min": Int % Range(15, 255, inclusive_end=True),
    "k_max": Int % Range(15, 255, inclusive_end=True),
    "k_step": Int % Range(2, 28, inclusive_end=True),
    "no_mercy": Bool,
    "bubble_level": Int % Range(0, 2, inclusive_end=True),
    "prune_level": Int % Range(0, 3, inclusive_end=True),
    "prune_depth": Int % Range(1, None),
    "disconnect_ratio": Float % Range(0, 1, inclusive_end=True),
    "low_local_ratio": Float % Range(0, 1, inclusive_end=True),
    "max_tip_len": Int % Range(1, None),
    "cleaning_rounds": Int % Range(1, None),
    "no_local": Bool,
    "kmin_1pass": Bool,
    "memory": Float % Range(0, None),
    "mem_flag": Int % Range(0, None),
    "num_cpu_threads": Int % Range(1, None),
    "no_hw_accel": Bool,
    "min_contig_len": Int,
}
# fmt: off
megahit_param_descriptions = {
    "presets": "Override a group of parameters. Possible values: 'meta-sensitive', "
               "'meta-large'.",
    "min_count": "Minimum multiplicity for filtering (k_min+1)-mers. Default: 2",
    "k_list": "List of kmer size - all must be odd with an increment <= 28. "
              "Default: [21,29,39,59,79,99,119,141]",
    "k_min": "Minimum kmer size (<= 255), must be odd number. "
             "Default: 21. Overrides k_list.",
    "k_max": "Maximum kmer size (<= 255), must be odd number. "
             "Default: 141. Overrides k_list.",
    "k_step": "Increment of kmer size of each iteration (<= 28), must be even number. "
              "Default: 12. Overrides k_list.",
    "no_mercy": "Do not add mercy kmers.",
    "bubble_level": "Intensity of bubble merging, 0 to disable. Default: 2.",
    "prune_level": "Strength of low depth pruning. Default: 2.",
    "prune_depth": "Remove unitigs with avg kmer depth less than this value. "
                   "Default: 2.",
    "disconnect_ratio": "Disconnect unitigs if its depth is less than this ratio times "
                        "the total depth of itself and its siblings. Default: 0.1.",
    "low_local_ratio": "Remove unitigs if its depth is less than this ratio times "
                       "the average depth of the neighborhoods. Default: 0.2.",
    "max_tip_len": "Remove tips less than this value.",
    "cleaning_rounds": "Number of rounds for graph cleanning. Default: 5.",
    "no_local": "Disable local assembly.",
    "kmin_1pass": "Use 1pass mode to build SdBG of k_min.",
    "memory": "Max memory in byte to be used in SdBG construction (if set between 0-1, "
              "fraction of the machine's total memory). Default: 0.9.",
    "mem_flag": "SdBG builder memory mode. 0: minimum; 1: moderate; "
    "others: use all memory specified by '-m/--memory'. Default: 1.",
    "num_cpu_threads": "Number of CPU threads. Default: # of logical processors.",
    "no_hw_accel": "Run MEGAHIT without BMI2 and POPCNT hardware instructions.",
    "min_contig_len": "Minimum length of contigs to output. Default: 200.",
}
# fmt: on
spades_params = {
    "isolate": Bool,
    "sc": Bool,
    "meta": Bool,
    "bio": Bool,
    "corona": Bool,
    "plasmid": Bool,
    "metaviral": Bool,
    "metaplasmid": Bool,
    "only_assembler": Bool,
    "careful": Bool,
    "disable_rr": Bool,
    "threads": Int % Range(1, None),
    "memory": Int % Range(1, None),
    "k": List[Int % Range(1, 128, inclusive_end=False)],
    "cov_cutoff": Float % Range(0, 1, inclusive_start=False)
    | Str % Choices(["auto", "off"]),
    "phred_offset": Int,
    "debug": Bool,
}
# fmt: off
spades_param_descriptions = {
    "isolate": "This flag is highly recommended for high-coverage isolate and "
               "multi-cell data.",
    "sc": "This flag is required for MDA (single-cell) data.",
    "meta": "This flag is required for metagenomic data.",
    "bio": "This flag is required for biosyntheticSPAdes mode.",
    "corona": "This flag is required for coronaSPAdes mode.",
    "plasmid": "Runs plasmidSPAdes pipeline for plasmid detection.",
    "metaviral": "Runs metaviralSPAdes pipeline for virus detection.",
    "metaplasmid": "Runs metaplasmidSPAdes pipeline for plasmid detection in "
                   "metagenomic datasets (equivalent for --meta --plasmid).",
    "only_assembler": "Runs only assembling (without read error correction).",
    "careful": "Tries to reduce number of mismatches and short indels.",
    "disable_rr": "Disables repeat resolution stage of assembling.",
    "threads": "Number of threads. Default: 16.",
    "memory": "RAM limit for SPAdes in Gb (terminates if exceeded). Default: 250.",
    "k": "List of k-mer sizes (must be odd and less than 128). Default: 'auto'.",
    "cov_cutoff": "Coverage cutoff value (a positive float number, or 'auto', or "
                  "'off'). Default: 'off'.",
    "phred_offset": "PHRED quality offset in the input reads (33 or 64). "
                    "Default: auto-detect.",
    "debug": "Runs SPAdes in debug mode.",
}
# fmt: on
quast_params = {
    # TODO: add eukaryote, fungal and large when alignment
    #  to reference is supported
    "min_contig": Int % Range(1, None),
    "threads": Int % Range(1, None),
    "k_mer_stats": Bool,
    "k_mer_size": Int % Range(1, None),
    "contig_thresholds": List[Int % Range(0, None)],
}
# fmt: off
quast_param_descriptions = {
    "min_contig": "Lower threshold for contig length. Default: 500.",
    "threads": "Maximum number of parallel jobs. Default: 25% of CPUs. "
               "Currently disabled - only 1 CPU is supported.",
    "k_mer_stats": "Compute k-mer-based quality metrics (recommended for large "
                   "genomes). This may significantly increase memory and time "
                   "consumption on large genomes.",
    "k_mer_size": "Size of k used in k-mer-stats. Default: 101.",
    "contig_thresholds": "List of contig length thresholds. "
                         "Default: 0,1000,5000,10000,25000,50000.",
}
# fmt: on
iss_params = {
    "sample_names": List[Str],
    "n_genomes": Int % Range(1, None),
    "ncbi": List[Str % Choices(["bacteria", "viruses", "archaea"])],
    "n_genomes_ncbi": List[Int % Range(1, None)],
    "abundance": Str
    % Choices(
        ["uniform", "halfnormal", "exponential", "lognormal", "zero-inflated-lognormal"]
    ),
    "coverage": Str
    % Choices(["halfnormal", "exponential", "lognormal", "zero-inflated-lognormal"]),
    "model": Str % Choices(["HiSeq", "NovaSeq", "MiSeq"]),
    "n_reads": Int % Range(1, None),
    "mode": Str % Choices(["kde", "basic", "perfect"]),
    "gc_bias": Bool,
    "cpus": Int % Range(1, None),
    "debug": Bool,
    "seed": Int % Range(1, None),
}
# fmt: off
iss_param_descriptions = {
    "sample_names": "List of sample names that should be generated. ",
    "n_genomes": "How many genomes will be used for the simulation. "
                 "Only required when genome sequences are provided. Default: 10.",
    "ncbi": "Download input genomes from NCBI. Can be bacteria, viruses, archaea "
            "or a combination of the three.",
    "n_genomes_ncbi": "How many genomes will be downloaded from NCBI. If more than "
                      "one kingdom is set with --ncbi, multiple values are necessary.",
    "abundance": "Abundance distribution. Default: lognormal.",
    "coverage": "Coverage distribution.",
    "model": "Error model. Use either of the precomputed models when --mode "
             "set to 'kde'. Default: HiSeq.",
    "n_reads": "Number of reads to generate. Default: 1000000.",
    "mode": "Error model. If not specified, using kernel density estimation. "
            "Default: kde.",
    "gc_bias": "If set, may fail to sequence reads with abnormal GC content.",
    "cpus": "Number of cpus to use. Default: 2.",
    "debug": " Enable debug logging.",
    "seed": "Seed for all the random number generators.",
}
# fmt: on
bowtie2_indexing_params = {
    "large_index": Bool,
    "debug": Bool,
    "sanitized": Bool,
    "verbose": Bool,
    "noauto": Bool,
    "packed": Bool,
    "bmax": Int % Range(1, None),
    "bmaxdivn": Int % Range(1, None),
    "dcv": Int % Range(1, None),
    "nodc": Bool,
    "offrate": Int % Range(0, None),
    "ftabchars": Int % Range(1, None),
    "threads": Int % Range(1, None),
    "seed": Int % Range(1, None),
}
# fmt: off
bowtie2_indexing_param_descriptions = {
    "large_index": "Force generated index to be 'large', even if ref has "
                   "fewer than 4 billion nucleotides.",
    "debug": "Use the debug binary; slower, assertions enabled.",
    "sanitized": "Use sanitized binary; slower, uses ASan and/or UBSan.",
    "verbose": "Log the issued command.",
    "noauto": " Disable automatic -p/--bmax/--dcv memory-fitting.",
    "packed": "Use packed strings internally; slower, less memory.",
    "bmax": "Max bucket sz for blockwise suffix-array builder.",
    "bmaxdivn": "Max bucket sz as divisor of ref len. Default: 4.",
    "dcv": "Diff-cover period for blockwise. Default: 1024.",
    "nodc": "Disable diff-cover (algorithm becomes quadratic).",
    "offrate": "SA is sampled every 2^<int> BWT chars. Default: 5.",
    "ftabchars": "# of chars consumed in initial lookup. Default: 10.",
    "threads": "# of CPUs.",
    "seed": "Seed for random number generator.",
}
# fmt: on
bowtie2_mapping_params = {
    "skip": Int % Range(1, None),
    "qupto": Int % Range(1, None),
    "trim5": Int % Range(0, None),
    "trim3": Int % Range(0, None),
    "trim_to": Str,
    "phred33": Bool,
    "phred64": Bool,
    "mode": Str % Choices(["local", "global"]),
    "sensitivity": Str % Choices(["very-fast", "fast", "sensitive", "very-sensitive"]),
    "n": Int % Range(0, 1, inclusive_end=True),
    "len": Int % Range(1, None),
    "i": Str,
    "n_ceil": Str,
    "dpad": Int % Range(0, None),
    "gbar": Int % Range(0, None),
    "ignore_quals": Bool,
    "nofw": Bool,
    "norc": Bool,
    "no_1mm_upfront": Bool,
    "end_to_end": Bool,
    "local": Bool,
    "ma": Int % Range(0, None),
    "mp": Str,
    "np": Int % Range(0, None),
    "rdg": Str,
    "rfg": Str,
    "k": Int % Range(0, None),
    "a": Bool,
    "d": Int % Range(0, None),
    "r": Int % Range(0, None),
    "minins": Int % Range(0, None),
    "maxins": Int % Range(1, None),
    "valid_mate_orientations": Str % Choices(["fr", "rf", "ff"]),
    "no_mixed": Bool,
    "no_discordant": Bool,
    "dovetail": Bool,
    "no_contain": Bool,
    "no_overlap": Bool,
    "offrate": Int % Range(0, None),
    "threads": Int % Range(1, None),
    "reorder": Bool,
    "mm": Bool,
    "seed": Int % Range(0, None),
    "non_deterministic": Bool,
}
# fmt: off
bowtie2_mapping_param_descriptions = {
    "skip": "Skip (i.e. do not align) the first <int> reads or pairs in the input.",
    "qupto": "Align the first <int> reads or read pairs from the input (after the "
             "-s/--skip reads or pairs have been skipped), then stop. "
             "Default: no limit.",
    "trim5": "Trim <int> bases from 5' (left) end of each read before alignment. "
             "Default: 0.",
    "trim3": "Trim <int> bases from 3' (right) end of each read before alignment "
             "Default: 0.",
    "trim_to": "Trim reads exceeding <int> bases. Bases will be trimmed from either "
               "the 3' (right) or 5' (left) end of the read. If the read end is not "
               "specified, bowtie2 will default to trimming from the 3' (right) end "
               "of the read. --trim-to and -trim3/-trim5 are mutually exclusive. "
               "The value of this parameter should have the following format: "
               "[3:|5:]<int>, e.g.: '5:120' if bases should be trimmed from 3' end or "
               "just '120' if the end is not specified.",
    "phred33": 'Input qualities are ASCII chars equal to the Phred quality plus 33. '
               'This is also called the "Phred+33" encoding, which is used by the '
               'very latest Illumina pipelines.',
    "phred64": 'Input qualities are ASCII chars equal to the Phred quality plus 64. '
               'This is also called the "Phred+64" encoding.',
    "mode": "bowtie2 alignment settings. See bowtie2 manual for more details.",
    "sensitivity": "bowtie2 alignment sensitivity. See bowtie2 manual for details.",
    "n": "Sets the number of mismatches to allowed in a seed alignment during "
         "multiseed alignment. Setting this higher makes alignment slower (often much "
         "slower) but increases sensitivity. Default: 0.",
    "len": "Sets the length of the seed substrings to align during multiseed "
           "alignment. Smaller values make alignment slower but more sensitive. "
           "Default: the --sensitive preset is used by default, which sets -L to "
           "22 and 20 in --end-to-end mode and in --local mode.",
    "i": 'Sets a function governing the interval between seed substrings to use during '
         'multiseed alignment. The value of this parameter should be provided as a '
         'comma-separated list, e.g.: "S,1,0.75". For details on how to set functions '
         'consult Bowtie 2 manual.',
    "n_ceil": 'Sets a function governing the maximum number of ambiguous characters '
              '(usually Ns and/or .s) allowed in a read as a function of read length. '
              'The value of this parameter should be provided as a comma-separated '
              'list, e.g.: "S,1,0.75". For details on how to set functions consult '
              'bowtie2 manual.',
    "dpad": '"Pads" dynamic programming problems by <int> columns on either side to '
            'allow gaps. Default: 15.',
    "gbar": "Disallow gaps within <int> positions of the beginning or end of the read. "
            "Default: 4.",
    "ignore_quals": "When calculating a mismatch penalty, always consider the quality "
                    "value at the mismatched position to be the highest possible, "
                    "regardless of the actual value. I.e. input is treated as though "
                    "all quality values are doesn't specify quality values (e.g. "
                    "in -f, -r, or -c modes).",
    "nofw": "If --nofw is specified, bowtie2 will not attempt to align unpaired reads "
            "to the forward (Watson) reference strand. In paired-end mode, pertains to "
            "the fragments. For more information, consult the Bowtie 2 manual.",
    "norc": "If --norc is specified, bowtie2 will not attempt to align unpaired reads "
            "against the reverse-complement (Crick) reference strand. In paired-end "
            "mode, pertains to the fragments. For more information, consult the "
            "bowtie2 manual.",
    "no_1mm_upfront": "By default, Bowtie 2 will attempt to find either an exact or a "
                      "1-mismatch end-to-end alignment for the read before trying the "
                      "multiseed heuristic. Such alignments can be found very quickly, "
                      "and many short read alignments have exact or near-exact "
                      "end-to-end alignments. However, this can lead to unexpected "
                      "alignments when the user also sets options governing the "
                      "multiseed heuristic, like -L and -N. For instance, if the user "
                      "specifies -N 0 and -L equal to the length of the read, the user "
                      "will be surprised to find 1-mismatch alignments reported. "
                      "This option prevents Bowtie 2 from searching for 1-mismatch "
                      "end-to-endalignments before using the multiseed heuristic, "
                      "which leads to the expected behavior when combined with options "
                      "such as -L and -N. This comes at the expense of speed.",
    "end_to_end": 'In this mode, Bowtie 2 requires that the entire read align from '
                  'one end to the other, without any trimming (or "soft clipping") of '
                  'characters from either end. The match bonus --ma always equals 0 in '
                  'this mode, so all alignment scores are less than or equal to 0, and '
                  'the greatest possible alignment score is 0. This is mutually '
                  'exclusive with --local. --end-to-end is the default mode.',
    "local": 'In this mode, bowtie2 does not require that the entire read align from '
             'one end to the other. Rather, some characters may be omitted '
             '("soft clipped") from the ends in order to achieve the greatest possible '
             'alignment score. The match bonus --ma is used in this mode, and the best '
             'possible alignment score is equal to the match bonus (--ma) times the '
             'length of the read. Specifying --local and one of the presets (e.g. '
             '--local --very-fast) is equivalent to specifying the local version of '
             'the preset (--very-fast-local). This is mutually exclusive with '
             '--end-to-end. --end-to-end is the default mode.',
    "ma": "Sets the match bonus. In --local mode <int> is added to the alignment score "
          "for each position where a read character aligns to a reference character "
          "and the characters match. Not used in --end-to-end mode. Default: 2.",
    "mp": "Sets the maximum (MX) and minimum (MN) mismatch penalties, both integers. "
          "A number less than or equal to MX and greater than or equal to MN is "
          "subtracted from the alignment score for each position where a read "
          "character aligns to a reference character, the characters do not match, "
          "and neither is an N. If --ignore-quals is specified, the number subtracted "
          "quals MX. Otherwise, the number subtracted is MN + floor((MX-MN)(MIN(Q,40.0"
          ")/40.0)) where Q is the Phred quality value. The value of this parameter "
          "should be provided as a comma-separated list of two integers. Default: 6,2.",
    "np": "Sets penalty for positions where the read, reference, or both, contain an "
          "ambiguous character such as N. Default: 1.",
    "rdg": "Sets the read gap open (<int1>) and extend (<int2>) penalties. A read gap "
           "of length N gets a penalty of <int1> + N * <int2>. The value of this "
           "parameter should be provided as a comma-separated list of two integers. "
           "Default: 5,3.",
    "rfg": "Sets the reference gap open (<int1>) and extend (<int2>) penalties. "
           "A reference gap of length N gets a penalty of <int1> + N * <int2>. The "
           "value of this parameter should be provided as a comma-separated list of "
           "two integers. Default: 5,3.",
    "k": "When -k is specified, however, bowtie2 searches for at most <int> distinct, "
         "valid alignments for each read. The search terminates when it can't find "
         "more distinct valid alignments, or when it finds <int>, whichever happens "
         "first. All alignments found are reported in descending order by alignment "
         "score. For more information, consult the bowtie2 manual.",
    "a": "Like -k but with no upper limit on number of alignments to search for. -a "
         "is mutually exclusive with -k.",
    "d": 'Up to <int> consecutive seed extension attempts can "fail" before bowtie2 '
         'moves on, using the alignments found so far. A seed extension "fails" if it '
         'does not yield a new best or a new second-best alignment. This limit is '
         'automatically adjusted up when -k or -a are specified. Default: 15.',
    "r": '<int> is the maximum number of times Bowtie 2 will "re-seed" reads with '
         'repetitive seeds. When "re-seeding," bowtie2 simply chooses a new set of '
         'reads (same length, same number of mismatches allowed) at different offsets '
         'and searches for more alignments. A read is considered to have repetitive '
         'seeds if the total number of seed hits divided by the number of seeds that '
         'aligned at least once is greater than 300. Default: 2.',
    "minins": "The minimum fragment length for valid paired-end alignments. "
              "Default: 0.",
    "maxins": "The maximum fragment length for valid paired-end alignments. "
              "Default: 500.",
    "valid_mate_orientations": "The upstream/downstream mate orientations for a valid "
                               "paired-end alignment against the forward reference "
                               "strand. For more details consult the bowtie2 manual.",
    "no_mixed": "By default, when bowtie2 cannot find a concordant or discordant "
                "alignment for a pair, it then tries to find alignments for the "
                "individual mates. This option disables that behavior.",
    "no_discordant": "By default, bowtie2 looks for discordant alignments if it "
                     "cannot find any concordant alignments. This option disables "
                     "that behavior.",
    "dovetail": 'If the mates "dovetail", that is if one mate alignment extends past '
                'the beginning of the other such that the wrong mate begins upstream, '
                'consider that to be concordant.',
    "no_contain": "If one mate alignment contains the other, consider that "
                  "to be non-concordant.",
    "no_overlap": "If one mate alignment overlaps the other at all, consider that "
                  "to be non-concordant.",
    "offrate": "Override the offrate of the index with <int>. If <int> is greater than "
               "the offrate used to build the index, then some row markings are "
               "discarded when the index is read into memory.",
    "threads": "Launch <int>> parallel search threads. Threads will run on separate "
               "processors/cores and synchronize when parsing reads and outputting "
               "alignments. Searching for alignments is highly parallel, and speedup "
               "is close to linear. Increasing -p increases Bowtie 2's memory "
               "footprint. Default: 1.",
    "reorder": "Guarantees that output SAM records are printed in an order "
               "corresponding to the order of the reads in the original input file, "
               "even when --threads is set greater than 1.",
    "mm": "Use memory-mapped I/O to load the index, rather than typical file I/O. "
          "Memory-mapping allows many concurrent bowtie processes on the same computer "
          "to share the same memory image of the index.",
    "seed": "Use <int> as the seed for pseudo-random number generator. Default: 0.",
    "non_deterministic": "If specified, Bowtie 2 re-initializes its pseudo-random "
                         "generator for each read using the current time. "
                         "Default: off.",
}
# fmt: on
