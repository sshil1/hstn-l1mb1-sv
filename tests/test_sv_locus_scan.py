import tempfile, os, gzip, sv_locus_scan as m
def w(name,txt,gz=False):
    p=os.path.join(tempfile.mkdtemp(),name)
    (gzip.open(p,'wt') if gz else open(p,'w')).write(txt); return p
def test_vcf_hit_and_size_flag():
    p=w('a.vcf.gz','##x\n#CHROM\tPOS\tID\tREF\tALT\tQ\tF\tINFO\n6\t85464318\t.\tG\t'+'A'*729+'\t.\t.\tSVTYPE=INS;SVLEN=728\n6\t10\t.\tG\tA\t.\t.\t.\n',True)
    r=m.scan([p],85459547,85467444,5000); assert len(r)==1 and r[0][6]=='SIZE_MATCH_600-850'
def test_vcf_without_svlen_uses_allele_length():
    p=w('b.vcf','6\t85464318\t.\tG\t'+'A'*729+'\t.\t.\tSVTYPE=INS\n'); r=m.scan([p],85459547,85467444,0); assert r and r[0][4]==728
def test_bed_outside_window_and_other_chrom():
    p=w('c.bed','6\t1000\t2000\n7\t85464318\t85464320\n6\t85470000\t85470100\n'); assert m.scan([p],85459547,85467444,1000)==[]
def test_bed_flank_hit():
    p=w('d.tsv','NC_037333.1\t85468000\t85468050\tDEL\n'); assert len(m.scan([p],85459547,85467444,1000))==1
def test_malformed_lines_do_not_crash():
    p=w('e.tsv','x\n6\tabc\tdef\n\n'); assert m.scan([p],1,2,0)==[]
