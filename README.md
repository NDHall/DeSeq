
##DeSeqTrinotateMerger.py  


DeSeqTrinotateMerger.py is written for python3, it mergers DeSeq csv output with Trinontate tsv output on the transcript_id field into a tsv file. It may be useful to know that all terms associated with a transcript_id are turned into 1 entry. For example, a transcript containing a gag-pol domain may have 4 separate entries assoicated with its transcript_id. DeSeqTrinotateMerger.py will report these as 1 entry.

##Usage:

    python3 DeSeqTrinotateMerger.py -d <DeSeq.csv> -t <Trinotate.tsv> -o <output_name.tsv>

    -d or --deseq     Comma delimited DeSeq.csv 
    -t or --trin      Tab delimited Trinontate file. This is the default for 
                      Trinontate and often has .xls extension.
    -o or --out       Full out file name. The out file will be tab delimited
                      this is to preserve the hiearchy of division employed 
                      within Trinotate's native output. 


     
