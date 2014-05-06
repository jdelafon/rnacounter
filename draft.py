import pysam
import os, sys
import itertools
from numpy import asarray


def parse_gtf(row):
    #fields = ['chr','source','name','start','end','score','strand','frame','attributes']
    def _score(x):
        if str(x) == '.': return 0
        else: return float(x)
    def _strand(x):
        smap = {'+':1, 1:1, '-':-1, -1:-1, '.':0, 0:0}
        return smap[x]
    if not row: return
    row = row.strip().split("\t")
    if len(row) < 9:
        raise ValueError("\"Attributes\" field required in GFF.")
    info = (row[0], int(row[3])-1, int(row[4]), row[2], _score(row[5]),_strand(row[6])) # bed fields
    attrs = (x.strip().split() for x in row[8].split(';'))  # {gene_id: "AAA", ...}
    attrs = dict((x[0],x[1].strip("\"")) for x in attrs)
    return info,attrs

class GenomicObject(object):
    def __init__(self, id='',gene_id='',gene_name='',chrom='',start=0,end=0,
                 name='',score=0.0,strand=0,length=0,seq=''):
        self.id = id
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand
        self.length = length
        self.seq = seq  # sequence
    def __and__(self,other):
        """The intersection of two GenomicObjects"""
        assert self.chrom==other.chrom, "Cannot add features from different chromosomes"
        selfid = (self.id,) if isinstance(self.id,int) else self.id
        otherid = (other.id,) if isinstance(other.id,int) else other.id
        return self.__class__(
            id = selfid + otherid,
            gene_id = '|'.join(set([self.gene_id, other.gene_id])),
            gene_name = '|'.join(set([self.gene_name, other.gene_name])),
            chrom = self.chrom,
            start = max(self.start, other.start),
            end = min(self.end, other.end),
            name = '|'.join(set([self.name, other.name])),
            #score = self.score + other.score,
            strand = (self.strand + other.strand)/2,
            #length = min(self.end, other.end) - max(self.start, other.start),
        )
    def __repr__(self):
        return "<%s (%d-%d) %s>" % (self.name,self.start,self.end,self.gene_name)

class Exon(GenomicObject):
    def __init__(self, transcripts=set(), **args):
        GenomicObject.__init__(self, **args)
        self.transcripts = transcripts   # list of transcripts it is contained in
        self.length = self.end - self.start
    def __and__(self,other):
        E = GenomicObject.__and__(self,other)
        E.transcripts = set(self.transcripts) | set(other.transcripts)
        return E

#class Gene(GenomicObject):
#    def __init__(self, exons=set(),transcripts=set(), **args):
#        GenomicObject.__init__(self, **args)
#        self.exons = exons               # list of exons contained
#        self.transcripts = transcripts   # list of transcripts contained

#class Transcript(GenomicObject):
#    def __init__(self, exons=[], **args):
#        GenomicObject.__init__(self, **args)
#        self.exons = exons               # list of exons it contains

def _inter(feats):
    """The intersection of a list of GenomicObjects"""
    return reduce(lambda x,y: x&y, feats)

def cobble(exons):
    ends = [(e.start,1,e) for e in exons] + [(e.end,0,e) for e in exons]
    ends.sort()
    active_exons = []
    cobbled = []
    for i in xrange(len(ends)-1):
        a = ends[i]
        b = ends[i+1]
        if a[1]==1:
            active_exons.append(a[2])
        elif a[1]==0:
            active_exons.remove(a[2])
        if len(active_exons)==0: continue
        if a[0]==b[0]: continue
        e = _inter(active_exons)
        e.start = a[0]; e.end = b[0];
        cobbled.append(e)
    return cobbled


######################################################################


def rnacount(bamname, annotname):
    sam = pysam.Samfile(bamname, "rb")
    annot = open(annotname, "r")

    row = annot.readline().strip()
    info,attrs = parse_gtf(row)
    chrom = info[0]
    start = int(info[1])-1
    lastchrom = chrom
    exonnr = 1

    while row:

        # Load all GTF exons of a chromosome in memory and sort
        chrexons = []
        while chrom == lastchrom:  # start <= lastend and
            info,attrs = parse_gtf(row)
            chrom = info[0]
            start = info[1]
            end = info[2]
            if (info[3] == "exon") and (end-start > 1):
                chrexons.append((info,attrs))
            row = annot.readline().strip()
            if not row:
                break
        lastchrom = chrom
        chrexons.sort(key=lambda x: (x[0][1],x[0][2]))  # sort w.r.t. start,end
        print ">> Chromosome", chrom

        lastend = sys.maxint
        lastgeneid = ''
        chunkexons = []
        for info,attrs in chrexons:
            start = info[1]
            end = info[2]
            gene_id = attrs['gene_id']
            # Store all overlapping intervals from a disjoint chunk
            if (start <= lastend) or (gene_id == lastgeneid):
            #if (start <= lastend):
                chunkexons.append((info,attrs))
                lastend = end
                lastgeneid = gene_id
                continue
            # Process the chunk of exons
            else:
                exons = []
                # Regroup info about the same exon into a single Exon instance
                for key,group in itertools.groupby(chunkexons, lambda x:x[0][1:3]):
                    # Group w.r.t. (start,end)
                    chunktrans = set()
                    for info,attrs in group:
                        chunktrans.add(attrs['transcript_id'])
                    # keep the last "attrs" that is common to all
                    # - unless an exon with the exact same coords is in two genes on diff strands...
                    exon_id = attrs.get('exon_id', 'E%d'%exonnr)
                    exons.append(Exon(id=exonnr, gene_id=attrs['gene_id'], gene_name=attrs['gene_name'],
                                      transcripts=chunktrans, chrom=info[0], start=info[1], end=info[2],
                                      name=exon_id, score=info[4], strand=info[5]))
                    exonnr += 1

                chunkexons = []
                lastend = sys.maxint

                # Cobble all these intervals
                pieces = cobble(exons)

                # Deduce the transcripts structure and filter too similar transcripts,
                # i.e. made of the same exons up to 200bp.
                t2e = {}                               # map {transcript: [exons ids]}
                for p in pieces:
                    if p.length < 200: continue        # filter out cobbled pieces of less that 200bp
                    for tx in p.transcripts:
                        t2e.setdefault(tx,[]).append(p.id)
                e2t = {}
                for t,e in t2e.iteritems():            # map {(exons ids combination): [transcripts]}
                    es = tuple(sorted(e))              # combination of exon indices
                    e2t.setdefault(es,[]).append(t)
                # replace too similar transcripts by the first of the list, arbitrarily
                transcripts = set()  # full list of remaining transcripts
                tx_replace = dict((badt,tlist[0]) for tlist in e2t.values() for badt in tlist[1:] if len(tlist)>1)
                for p in pieces:
                    filtered = set([tx_replace.get(t,t) for t in p.transcripts])
                    transcripts |= filtered
                    p.transcripts = list(filtered)
                transcripts = list(transcripts)
                # remake the transcript-exon mapping
                t2e = {}
                for p in pieces:
                    for tx in p.transcripts:
                        t2e.setdefault(tx,[]).append(p.name)

                # Build the matrix :: lines are exons, columns are transcripts,
                # so that A[i,j]=1 means "transcript j contains exon i".
                Avals = asarray([[int(t in e.transcripts) for t in transcripts] for e in pieces])

                if exons[0].gene_name == "Gapdh":
                    print t2e
                    print exons
                    print transcripts
                    print Avals
                #if len(Avals)>0:
                #    print Avals

                #break    # 1 feat per chromosome
        #break    # 1 chromosome

    annot.close()
    sam.close


######################################################################

bamname = "testfiles/test_input_EV_chr1p_flash_htsstation.bam"
annotname = "testfiles/mm9.gtf"

rnacount(bamname,annotname)


