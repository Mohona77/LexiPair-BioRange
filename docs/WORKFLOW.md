# WORKFLOW (Slide 11 METHOD 1-7 -> code)
1 ref (tools/generate_datasets.py -> data/reference_*) / 2 haps (same -> data/primary_haplotypes.fasta) / 3 encode (include/lexipair/kmer_codec.h) / 4 sort+run (lexipair_index.h) / 5 cache layout (kmer_entry.h 16B, 4/line) / 6a static (apps/static_query) / 6b LSM (lsm_index.h) / 6c FC (fractional_cascade.h) / 7 validate (tests/ + results/validation.txt).
