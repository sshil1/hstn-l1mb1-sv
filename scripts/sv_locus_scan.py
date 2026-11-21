#!/usr/bin/env python3
"""Scan SV supplement files for records near the HSTN intron-1 locus (BTA6).
Read-only. Accepts .vcf/.vcf.gz, .bed/.tsv/.txt/.csv and .xlsx (needs openpyxl).
Usage: sv_locus_scan.py [--start N --end N --flank N] FILE [FILE ...]
Output: TSV to stdout; one row per hit, with a note on size match (600-850 bp)."""
import sys, gzip, csv, re, argparse, io
CHR6={'6','chr6','Chr6','NC_037333.1'}  # add other chr6 aliases used by a given supplement here
def opener(p):
    return gzip.open(p,'rt',errors='replace') if p.endswith('.gz') else open(p,'r',errors='replace')
def rows_from_vcf(f):
    for ln in f:
        if ln.startswith('#') or not ln.strip(): continue
        c=ln.rstrip('\n').split('\t')
        if len(c)<5: continue
        info=c[7] if len(c)>7 else ''
        m=re.search(r'(?:^|;)END=(\d+)',info); s=re.search(r'(?:^|;)SVLEN=(-?\d+)',info); t=re.search(r'(?:^|;)SVTYPE=([A-Za-z]+)',info)
        pos=int(c[1]); end=int(m.group(1)) if m else pos+max(len(c[3]),1)-1
        size=abs(int(s.group(1))) if s else abs(len(c[4])-len(c[3]))
        yield c[0],pos,end,size,(t.group(1) if t else ''),ln.strip()[:160]
def rows_from_table(f,delim=None):
    text=f.read(); delim=delim or ('\t' if text.count('\t')>=text.count(',') else ',')
    for c in csv.reader(io.StringIO(text),delimiter=delim):
        if len(c)<3: continue
        chrom=c[0].strip()
        nums=[x for x in c[1:] if re.fullmatch(r'\d+',x.strip())]
        if chrom not in CHR6 or len(nums)<2: continue
        s,e=int(nums[0]),int(nums[1])
        txt='\t'.join(c)[:160]
        sz=re.search(r'(?:INS|DEL)[^\d\-]*(\d{2,5})',txt)
        yield chrom,s,e,int(sz.group(1)) if sz else abs(e-s),'',txt
def rows_from_xlsx(p):
    import openpyxl
    wb=openpyxl.load_workbook(p,read_only=True,data_only=True)
    for ws in wb:
        for r in ws.iter_rows(values_only=True):
            r=['' if v is None else str(v) for v in r]
            if len(r)<3: continue
            if r[0].strip() in CHR6:
                nums=[x for x in r[1:] if re.fullmatch(r'\d+(\.0)?',x.strip())]
                if len(nums)>=2: yield r[0],int(float(nums[0])),int(float(nums[1])),abs(int(float(nums[1]))-int(float(nums[0]))),ws.title,'\t'.join(r)[:160]
def scan(paths,start,end,flank):
    out=[]
    for p in paths:
        if p.lower().endswith('.xlsx'): it=rows_from_xlsx(p)
        else:
            f=opener(p); it=rows_from_vcf(f) if re.search(r'\.vcf(\.gz)?$',p.lower()) else rows_from_table(f)
        for chrom,s,e,size,t,raw in it:
            if chrom in CHR6 and e>=start-flank and s<=end+flank:
                out.append((p,chrom,s,e,size,t,'SIZE_MATCH_600-850' if 600<=size<=850 else '',raw))
    return out

BIG=re.compile(r'(?<![\d.])(\d{6,10})(?:\.0)?(?![\d])')
IDCHR=re.compile(r'(?:^|[^\w.])(?:chr)?6[-:_]\d{6,}')
def all_rows(p):
    """Yield (sheet, [cells as str]) for any supported file; layout-agnostic."""
    if p.lower().endswith('.xlsx'):
        import openpyxl
        wb=openpyxl.load_workbook(p,read_only=True,data_only=True)
        for ws in wb:
            for r in ws.iter_rows(values_only=True):
                yield ws.title,['' if v is None else str(v) for v in r]
    else:
        with opener(p) as f:
            for ln in f:
                if ln.startswith('#') or not ln.strip(): continue
                yield '',(ln.rstrip('\n').split('\t') if '\t' in ln else next(csv.reader([ln.rstrip('\n')]),[]))
def scan_loose(paths,start,end,flank):
    """Layout-agnostic scan: flag any row with an 6-10 digit integer inside the window.
    LOOSE_CHR6 = row also carries a chr6 label (cell or 6-POS/chr6:POS id) or two coordinates spanning the window;
    LOOSE_COORD_ONLY = coordinate in window but no chr6 label (needs manual check; may be another chromosome)."""
    lo,hi=start-flank,end+flank; out=[]
    for p in paths:
        try:
            for sheet,cells in all_rows(p):
                text='\t'.join(cells)
                # coordinates come from purely numeric cells (so gene ids like LOC100336682 are ignored) or chr6-POS ids
                nums=[int(float(c.strip())) for c in cells if re.fullmatch(r'\d{6,10}(?:\.0)?',c.strip())]
                nums+= [int(x) for x in re.findall(r'(?:^|[^\w.])(?:chr)?6[-:_](\d{6,10})',text)]
                idm=bool(IDCHR.search(text))
                chr6=any(c.strip() in CHR6 for c in cells) or idm or 'NC_037333.1' in text
                inwin=[n for n in nums if lo<=n<=hi]
                span=chr6 and len(nums)>=2 and min(nums)<=hi and max(nums)>=lo
                if inwin or span:
                    out.append((p,'LOOSE_CHR6' if chr6 else 'LOOSE_COORD_ONLY',sheet,text[:200]))
        except Exception as ex:
            out.append((p,'ERROR',type(ex).__name__,str(ex)[:200]))
    return out
def report(paths):
    """Per-file parse diagnostics so that '0 hits' can be judged: rows, sheets, strict chr6 rows, error."""
    rep=[]
    for p in paths:
        d={'file':p,'rows':0,'sheets':[],'chr6_rows':0,'error':''}
        try:
            for sheet,cells in all_rows(p):
                d['rows']+=1
                if sheet and sheet not in d['sheets']: d['sheets'].append(sheet)
                if any(c.strip() in CHR6 for c in cells) or IDCHR.search('\t'.join(cells)): d['chr6_rows']+=1
        except Exception as ex:
            d['error']=type(ex).__name__+': '+str(ex)[:120]
        rep.append(d)
    return rep
def main(argv=None):
    ap=argparse.ArgumentParser(); ap.add_argument('--start',type=int,default=85459547); ap.add_argument('--end',type=int,default=85467444)
    ap.add_argument('--flank',type=int,default=5000); ap.add_argument('--loose',action='store_true',help='layout-agnostic scan (any column)'); ap.add_argument('--report',action='store_true',help='print per-file parse diagnostics to stderr'); ap.add_argument('files',nargs='+'); a=ap.parse_args(argv)
    if a.report:
        for d in report(a.files): print('REPORT\t%(file)s\trows=%(rows)d\tchr6_rows=%(chr6_rows)d\tsheets=%(sheets)s\t%(error)s'%d,file=sys.stderr)
    if a.loose:
        print('file\tflag\tsheet\trecord'); L=scan_loose(a.files,a.start,a.end,a.flank)
        for r in L: print('\t'.join(map(str,r)))
        print(f'# {len(L)} loose hit(s) in {a.start-a.flank}-{a.end+a.flank}; LOOSE_COORD_ONLY rows need manual chromosome check',file=sys.stderr); return
    print('file\tchrom\tstart\tend\tsize\ttype\tflag\trecord')
    n=0
    for r in scan(a.files,a.start,a.end,a.flank): print('\t'.join(map(str,r))); n+=1
    print(f'# {n} hit(s) within HSTN intron 1 +/- {a.flank} bp ({a.start}-{a.end}); coordinates are compared as written, build differences are not lifted over',file=sys.stderr)
if __name__=='__main__': main()
