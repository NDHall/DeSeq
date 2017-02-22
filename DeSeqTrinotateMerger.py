 
import sys
import getopt

class trinontate_entry:

    def __init__(self, gene_id, transcript_id, sprot_top_blastx_hit, trembl_top_blastx_hit, rnammer, prot_id, prot_coords, sprot_top_blastp_hit, trembl_top_blastp_hit, pfam, signalp,tmhh, eggnog, gene_ontology_blast, gene_ontology_pfam, transcript, peptide ):
        self.gene_id=str(gene_id)
        self.transcript_id=str(transcript_id)
        self.sprot_top_blastx_hit=sprot_top_blastx_hit
        self.trembl_top_blastx_hit=trembl_top_blastx_hit
        self.rnammer=rnammer
        self.prot_id=prot_id
        self.prot_coords=prot_coords
        self.sprot_top_blastp_hit=sprot_top_blastp_hit
        self.trembl_top_blastp_hit=trembl_top_blastp_hit
        self.pfam=pfam
        self.signalp=signalp
        self.tmhh=tmhh
        self.eggnog=eggnog
        self.gene_ontology_blast=gene_ontology_blast
        self.gene_ontology_pfam=gene_ontology_pfam
        self.transcript=transcript
        self.peptide=peptide
    def add_to_gene_id(self, atr,delim):
        if atr != self.gene_id or self.gene_id=="." or atr == ".":
            self.gene_id=str(self.gene_id)+delim+atr
    def add_to_transcript_id(self, atr,delim):
        if atr != self.transcript_id or self.transcript_id=="." or atr == ".":
            self.transcript_id=str(self.transcript_id)+delim+atr
    def add_to_sprot_top_blastx_hit(self, atr,delim):
        if atr != self.sprot_top_blastx_hit or self.sprot_top_blastx_hit=="." or atr == ".":
            self.sprot_top_blastx_hit=str(self.sprot_top_blastx_hit)+delim+atr
    def add_to_trembl_top_blastx_hit(self, atr,delim):
        if atr != self.trembl_top_blastx_hit or self.trembl_top_blastx_hit=="." or atr == ".":
            self.trembl_top_blastx_hit=str(self.trembl_top_blastx_hit)+delim+atr
    def add_to_rnammer(self, atr,delim):
        if atr != self.rnammer or self.rnammer=="." or atr == ".":
            self.rnammer=str(self.rnammer)+delim+atr
    def add_to_prot_id(self, atr,delim):
        if atr != self.prot_id or self.prot_id=="." or atr == ".":
            self.prot_id=str(self.prot_id)+delim+atr
    def add_to_prot_coords(self, atr,delim):
        if atr != self.prot_coords or self.prot_coords=="." or atr == ".":
            self.prot_coords=str(self.prot_coords)+delim+atr
    def add_to_sprot_top_blastp_hit(self, atr,delim):
        if atr != self.sprot_top_blastp_hit or self.sprot_top_blastp_hit=="." or atr == ".":
            self.sprot_top_blastp_hit=str(self.sprot_top_blastp_hit)+delim+atr
    def add_to_trembl_top_blastp_hit(self, atr,delim):
        if atr != self.trembl_top_blastp_hit or self.trembl_top_blastp_hit=="." or atr == ".":
            self.trembl_top_blastp_hit=str(self.trembl_top_blastp_hit)+delim+atr
    def add_to_pfam(self, atr,delim):
        if atr != self.pfam or self.pfam=="." or atr == ".":
            self.pfam=str(self.pfam)+delim+atr
    def add_to_signalp(self, atr,delim):
        if atr != self.signalp or self.signalp=="." or atr == ".":
            self.signalp=str(self.signalp)+delim+atr
    def add_to_tmhh(self, atr,delim):
        if atr != self.tmhh or self.tmhh=="." or atr == ".":
            self.tmhh=str(self.tmhh)+delim+atr
    def add_to_eggnog(self, atr,delim):
        if atr != self.eggnog or self.eggnog=="." or atr == ".":
            self.eggnog=str(self.eggnog)+delim+atr
    def add_to_gene_ontology_blast(self, atr,delim):
        if atr != self.gene_ontology_blast or self.gene_ontology_blast=="." or atr == ".":
            self.gene_ontology_blast=str(self.gene_ontology_blast)+delim+atr
    def add_to_gene_ontology_pfam(self, atr,delim):
        if atr != self.gene_ontology_pfam or self.gene_ontology_pfam=="." or atr == ".":
            self.gene_ontology_pfam=str(self.gene_ontology_pfam)+delim+atr
    def add_to_transcript(self, atr,delim):
        if atr != self.transcript or self.transcript=="." or atr == ".":
            self.transcript=str(self.transcript)+delim+atr
    def add_to_peptide(self, atr,delim):
        if atr != self.peptide or self.peptide=="." or atr == ".":
            self.peptide=str(self.peptide)+delim+atr

    def trinotate_entry_pile(self, gene_id, transcript_id, sprot_top_blastx_hit, trembl_top_blastx_hit, rnammer, prot_id, prot_coords, sprot_top_blastp_hit, trembl_top_blastp_hit, pfam, signalp,tmhh, eggnog, gene_ontology_blast, gene_ontology_pfam, transcript, peptide ):
        self.add_to_gene_id(gene_id , "`")
        self.add_to_transcript_id(transcript_id , "`")
        self.add_to_sprot_top_blastx_hit(sprot_top_blastx_hit , "`")
        self.add_to_trembl_top_blastx_hit(trembl_top_blastx_hit , "`")
        self.add_to_rnammer(rnammer , "`")
        self.add_to_prot_id(prot_id , "`")
        self.add_to_prot_coords(prot_coords , "`")
        self.add_to_sprot_top_blastp_hit(sprot_top_blastp_hit , "`")
        self.add_to_trembl_top_blastp_hit(trembl_top_blastp_hit , "`")
        self.add_to_pfam(pfam , "`")
        self.add_to_signalp(signalp , "`")
        self.add_to_tmhh(tmhh , "`")
        self.add_to_eggnog(eggnog , "`")
        self.add_to_gene_ontology_blast(gene_ontology_blast , "`")
        self.add_to_gene_ontology_pfam(gene_ontology_pfam , "`")
        self.add_to_transcript(transcript , "`")
        self.add_to_peptide(peptide , "`")


class deseq_out:
    def __init__(self,transcript_id, base_mean, log2_fc, lfc_se, stat, pvalue, padj ):
        self.transcript_id=str(transcript_id)
        self.base_mean=base_mean
        self.log2_fc=log2_fc
        self.lfc_se=lfc_se
        self.stat=stat
        self.pvalue=pvalue
        self.padj=padj




def clean_for_append(out):
    f=open(out,'w')
    f.write("")
    f.close()
    f=open(out,'a')
    return f

def main(argv):
    argv_dict={"out":"", "trin":"", "deseq":""}
    try:
        opts, args = getopt.getopt(argv,"h:t:d:o:",["trinnotate=","deseq=","out="])
    except getopt.error:
        print ("getopt failed")
        print ("\n\n\tDeSeqTrinontateMerger.py -t <trinnontate file> -d <deseq file> -o <out>\n\n")
        sys.exit(2)
    for opt , arg in opts:
        if opt == '-h':
            print("\n\n\tDeSeqTrinontateMerger.py -t <trinnontate file> -d <deseq file> -o <out>\n\n")
            sys.exit(0)
        elif opt in ("-d", "--deseq" ):
            argv_dict["deseq"]=arg
        elif opt in ("-t","--trinnontate"):
            argv_dict["trin"]=arg
        elif opt in ("-o","--out"):
            argv_dict["out"]=arg
    return argv_dict
                        


#==================================================

#        Function parses tab delimited trinontate
#        file and returns a dictionary with 
#        transcript_id column 2 as the the key
#        for each entry. 

#==================================================



def trinotate_parse(file):
    counter=0
    trin_dict={}
    with open(file) as f:
        content = f.readlines()
        for line in content:
            line=line.strip()
            line=line.split("\t")
            if counter == 0:
                assert("gene_id" in line[0] and  "transcript_id" in line[1] and  "sprot_Top_BLASTX_hit" in line[2] and  "TrEMBL_Top_BLASTX_hit" in line[3] and  "RNAMMER" in line[4] and  "prot_id" in line[5] and  "prot_coords" in line[6] and  "sprot_Top_BLASTP_hit" in line[7] and  "TrEMBL_Top_BLASTP_hit" in line[8] and  "Pfam" in line[9] and  "SignalP" in line[10] and  "TmHMM" in line[11] and  "eggnog" in line[12] and  "gene_ontology_blast" in line[13] and  "gene_ontology_pfam" in line[14] and  "transcript" in line[15] and  "peptide" in line[16]), "\n\n\t %r is not a tab delimited trinnontate file. Please make sure that the headings are not modified from the orginal output and that tab characters separate each column)" % file
                counter+=1
            else :

                if line[1] not in trin_dict:
                    trin_line = trinontate_entry(line[0], line[1], line[2], 
                                           line[3], line[4], line[5], 
                                           line[6], line[7], line[8], 
                                           line[9], line[10],line[11], 
                                           line[12],line[13],line[14], 
                                           line[15],line[16])
                                 
                    trin_dict[trin_line.transcript_id]=trin_line
                else:
                    trin_dict[line[1]].trinotate_entry_pile(line[0], line[1], line[2], 
                                           line[3], line[4], line[5], 
                                           line[6], line[7], line[8], 
                                           line[9], line[10],line[11], 
                                           line[12],line[13],line[14], 
                                           line[15],line[16])
                   
                
                counter+=1
    f.close() 
    return trin_dict

def deseq_parse_join(file, trin_dict, out):
    counter=0
    out_file=clean_for_append(out)
    out_file.write("\t".join(["transcipt_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "gene_id", "transcript_id", "sprot_Top_BLASTX_hit", "TrEMBL_Top_BLASTX_hit", "RNAMMER", "prot_id", "prot_coords", "sprot_Top_BLASTP_hit", "TrEMBL_Top_BLASTP_hit", "Pfam", "SignalP", "TmHMM", "eggnog", "gene_ontology_blast", "gene_ontology_pfam", "transcript", "peptide"])+"\n")

    with open(file) as f:
        content = f.readlines()
        for line in content:
            line=line.strip()
            line=line.split(",")
            assert(len(line)==7), "\n\n\t %r does not appear to be a csv file please check and make sure that you have the correct input file.\n\n." % file
            if counter == 0:
                assert("baseMean" in line[1] and "log2FoldChange" in line[2] and "lfcSE" in line[3] and "stat" in line[4] and "pvalue" in line[5] and "padj" in line[6]), "\n\n\y%r does not appear an unaltered DeSeq output. Please double check your file, and make sure that it is in csv form and that you have the correct file.\n\n" % file 
                counter+= 1 
            else: 
                deseq_line=deseq_out(line[0], line[1], line[2], line[3], line[4], line[5], line[6])
                assert(deseq_line.transcript_id in trin_dict), "\n\n\t%r does not have a matching entry in the trinnotate output.\n\n" % deseq_line.transcript_id
                trin_line=trin_dict[deseq_line.transcript_id]
                out_file.write("\t".join([deseq_line.transcript_id ,deseq_line.base_mean ,deseq_line.log2_fc ,deseq_line.lfc_se ,deseq_line.stat ,deseq_line.pvalue ,deseq_line.padj ,trin_line.gene_id, trin_line.transcript_id, trin_line.sprot_top_blastx_hit, trin_line.trembl_top_blastx_hit, trin_line.rnammer, trin_line.prot_id, trin_line.prot_coords, trin_line.sprot_top_blastp_hit, trin_line.trembl_top_blastp_hit, trin_line.pfam, trin_line.signalp, trin_line.tmhh, trin_line.eggnog, trin_line.gene_ontology_blast, trin_line.gene_ontology_pfam, trin_line.transcript, trin_line.peptide])+"\n")
    out_file.close()

if __name__== "__main__" :
    argv=main(sys.argv[1:])
    trin_dict=trinotate_parse(argv["trin"])
    deseq_parse_join(argv["deseq"], trin_dict, argv["out"])
    
    











            
    
