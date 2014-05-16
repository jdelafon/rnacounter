
import pysam
import os, sys
import itertools
from numpy import asarray, zeros
import numpy
import copy


Ecounter = itertools.count(1)  # to give unique ids to undefined exons, see parse_gtf()

def parse_gtf(row):
    # GTF fields = ['chr','source','name','start','end','score','strand','frame','attributes']
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
    if row[2] != 'exon':
        return False
    attrs = (x.strip().split() for x in row[8].split(';'))  # {gene_id: "AAA", ...}
    attrs = dict((x[0],x[1].strip("\"")) for x in attrs)
    exon_id = attrs.get('exon_id', 'E%d'%Ecounter.next())
    return Exon(id=exon_id, gene_id=attrs['gene_id'], gene_name=attrs['gene_name'],
                chrom=row[0], start=int(row[3])-1, end=int(row[4]),
                name=exon_id, score=_score(row[5]), strand=_strand(row[6]),
                transcripts=[attrs['transcript_id']])


def lsqnonneg(C, d, x0=None, tol=None, itmax_factor=3):
    """Linear least squares with nonnegativity constraints (NNLS), based on MATLAB's lsqnonneg function.

    ``(x,resnorm,res) = lsqnonneg(C,d)`` returns

    * the vector *x* that minimizes norm(d-Cx) subject to x >= 0
    * the norm of residuals *resnorm* = norm(d-Cx)^2
    * the residuals *res* = d-Cx

    :param x0: Initial point for x.
    :param tol: Tolerance to determine what is considered as close enough to zero.
    :param itmax_factor: Maximum number of iterations.

    :type C: nxm numpy array
    :type d: nx1 numpy array
    :type x0: mx1 numpy array
    :type tol: float
    :type itmax_factor: int
    :rtype: *x*: numpy array, *resnorm*: float, *res*: numpy array

    Reference: C.L. Lawson and R.J. Hanson, Solving Least-Squares Problems, Prentice-Hall, Chapter 23, p. 161, 1974.
    `<http://diffusion-mri.googlecode.com/svn/trunk/Python/lsqnonneg.py>`_
    """
    def norm1(x):
        return abs(x).sum().max()

    def msize(x, dim):
        s = x.shape
        if dim >= len(s): return 1
        else: return s[dim]

    def positive(x):
        """Set to zero all negative components of an array."""
        for i in range(len(x)):
            if x[i] < 0 : x[i] = 0
        return x

    eps = sys.float_info.epsilon
    if tol is None: tol = 10*eps*norm1(C)*(max(C.shape)+1)
    C = asarray(C)
    (m,n) = C.shape
    P = numpy.zeros(n)                       # set P of indices where ultimately x_j > 0 & w_j = 0
    Z = ZZ = numpy.arange(1, n+1)            # set Z of indices where ultimately x_j = 0 & w_j <= 0
    if x0 is None or any(x0 < 0): x = P
    else: x = x0
    resid = d - numpy.dot(C, x)
    w = numpy.dot(C.T, resid)                # n-vector C'(d-Cx), "dual" of x, gradient of (1/2)*||d-Cx||^2
    outeriter = it = 0
    itmax = itmax_factor*n
    # Outer loop "A" to hold positive coefficients
    while numpy.any(Z) and numpy.any(w[ZZ-1] > tol): # if Z is empty or w_j<0 for all j, terminate.
        outeriter += 1
        t = w[ZZ-1].argmax()                 # find index t s.t. w_t = max(w), w_t in Z. So w_t > 0.
        t = ZZ[t]
        P[t-1]=t                             # move the index t from set Z to set P
        Z[t-1]=0                             # Z becomes [0] if n=1
        PP = numpy.where(P != 0)[0]+1        # non-zero elements of P for indexing, +1 (-1 later)
        ZZ = numpy.where(Z != 0)[0]+1        # non-zero elements of Z for indexing, +1 (-1 later)
        CP = numpy.zeros(C.shape)
        CP[:, PP-1] = C[:, PP-1]             # CP[:,j] is C[:,j] if j in P, or 0 if j in Z
        CP[:, ZZ-1] = numpy.zeros((m, msize(ZZ, 1)))
        z=numpy.dot(numpy.linalg.pinv(CP), d)                 # n-vector solution of least-squares min||d-CPx||
        if isinstance(ZZ,numpy.ndarray) and len(ZZ) == 0:     # if Z = [0], ZZ = [] and makes it fail
            return (positive(z), sum(resid*resid), resid)
        else:
            z[ZZ-1] = numpy.zeros((msize(ZZ,1), msize(ZZ,0))) # define z_j := 0 for j in Z
        # Inner loop "B" to remove negative elements from z if necessary
        while numpy.any(z[PP-1] <= tol): # if z_j>0 for all j, set x=z and return to outer loop
            it += 1
            if it > itmax:
                max_error = z[PP-1].max()
                raise Exception('Exiting: Iteration count (=%d) exceeded\n \
                      Try raising the tolerance tol. (max_error=%d)' % (it, max_error))
            QQ = numpy.where((z <= tol) & (P != 0))[0]        # indices j in P s.t. z_j < 0
            alpha = min(x[QQ]/(x[QQ] - z[QQ]))                # step chosen as large as possible s.t. x remains >= 0
            x = x + alpha*(z-x)                               # move x by this step
            ij = numpy.where((abs(x) < tol) & (P <> 0))[0]+1  # indices j in P for which x_j = 0
            Z[ij-1] = ij                                      # Add to Z, remove from P
            P[ij-1] = numpy.zeros(max(ij.shape))
            PP = numpy.where(P != 0)[0]+1
            ZZ = numpy.where(Z != 0)[0]+1
            CP[:, PP-1] = C[:, PP-1]
            CP[:, ZZ-1] = numpy.zeros((m, msize(ZZ, 1)))
            z=numpy.dot(numpy.linalg.pinv(CP), d)
            z[ZZ-1] = numpy.zeros((msize(ZZ,1), msize(ZZ,0)))
        x = z
        resid = d - numpy.dot(C, x)
        w = numpy.dot(C.T, resid)
    return (x, sum(resid*resid), resid)


class GenomicObject(object):
    def __init__(self, id='',gene_id='',gene_name='',chrom='',start=0,end=0,
                 name='',score=0.0,count=0,rpk=0.0,strand=0,length=0,seq='',multiplicity=1):
        self.id = id
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        #self.score = score
        self.count = count
        self.rpk = rpk
        self.strand = strand
        self.length = length
        #self.seq = seq  # sequence
        self.multiplicity = multiplicity
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
            #start = max(self.start, other.start),
            #end = min(self.end, other.end),
            ##   name = '|'.join(set([self.name, other.name])),
            name = '|'.join([self.name, other.name]),
            #score = self.score + other.score,
            strand = (self.strand + other.strand)/2,
            #length = min(self.end, other.end) - max(self.start, other.start),
            multiplicity = self.multiplicity + other.multiplicity
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

class Transcript(GenomicObject):
    def __init__(self, exons=[], **args):
        GenomicObject.__init__(self, **args)
        self.exons = exons               # list of exons it contains

class Gene(GenomicObject):
    def __init__(self, exons=set(),transcripts=set(), **args):
        GenomicObject.__init__(self, **args)
        self.exons = exons               # list of exons contained
        self.transcripts = transcripts   # list of transcripts contained


def intersect_exons_list(feats, multiple=False):
    """The intersection of a list *feats* of GenomicObjects.
    If *multiple* is True, permits multiplicity: if the same exon E1 is
    given twice, there will be "E1|E1" parts. Otherwise pieces are unique."""
    if multiple is False:
        feats = list(set(feats))
    if len(feats) == 1:
        return copy.deepcopy(feats[0])
    else:
        return reduce(lambda x,y: x&y, feats)


def cobble(exons, multiple=False):
    """Split exons into non-overlapping parts.
    :param multiple: see intersect_exons_list()."""
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
        if len(active_exons)==0:
            continue
        if a[0]==b[0]:
            continue
        e = intersect_exons_list(active_exons)
        e.start = a[0]; e.end = b[0]; e.length = b[0]-a[0]
        cobbled.append(e)
    return cobbled


class Counter(object):
    def __init__(self, stranded=False):
        self.n = 0 # read count
        self.n_raw = 0 # read count, no NH flag
        self.n_ws = 0 # read count, wrong strand
        self.strand = 0 # exon strand
        if stranded:
            self.count_fct = self.count_stranded
        else:
            self.count_fct = self.count

    def __call__(self, alignment):
        self.count_fct(alignment)

    def count(self, alignment):
        NH = [1.0/t[1] for t in alignment.tags if t[0]=='NH']+[1]
        self.n += NH[0]
        self.n_raw += 1

    def count_stranded(self, alignment):
        NH = [1.0/t[1] for t in alignment.tags if t[0]=='NH']+[1]
        if self.strand == 1 and alignment.is_reverse == False \
        or self.strand == -1 and alignment.is_reverse == True:
            self.n += NH[0]
        else:
            self.n_ws += NH[0]


# Gapdh id: ENSMUSG00000057666
# Gapdh transcripts: ENSMUST00000147954, ENSMUST00000147954, ENSMUST00000118875
#                    ENSMUST00000073605, ENSMUST00000144205, ENSMUST00000144588


######################################################################

def process_chunk(ckexons, sam, chrom, lastend):
    """Distribute counts across transcripts and genes of a chunk *ckexons*
    of non-overlapping exons."""

    def toRPK(count,length):
        return 1000.0 * count / length
    def fromRPK(rpk,length):
        return length * rpk / 1000.

    #if ckexons[0].gene_name != "Gapdh": return 1


    #--- Convert chromosome name
    if chrom[:3] == "NC_" : pass
    try: int(chrom); chrom = "chr"+chrom
    except: pass


    #--- Regroup occurrences of the same Exon from a different transcript
    exons = []
    for key,group in itertools.groupby(ckexons, lambda x:x.id):
        # ckexons are sorted because chrexons were sorted by chrom,start,end
        exon0 = group.next()
        for g in group:
            exon0.transcripts.append(g.transcripts[0])
        exons.append(exon0)
    del ckexons


    #--- Cobble all these intervals
    pieces = cobble(exons)


    #--- Filter out too similar transcripts,
    # e.g. made of the same exons up to 100bp.
    t2e = {}                               # map {transcript: [pieces IDs]}
    for p in pieces:
        if p.length < 100: continue        # filter out cobbled pieces of less that read length
        for t in p.transcripts:
            t2e.setdefault(t,[]).append(p.id)
    e2t = {}
    for t,e in t2e.iteritems():
        es = tuple(sorted(e))              # combination of pieces indices
        e2t.setdefault(es,[]).append(t)    # {(pieces IDs combination): [transcripts with same struct]}
    # Replace too similar transcripts by the first of the list, arbitrarily
    transcript_ids = set()  # full list of remaining transcripts
    tx_replace = dict((badt,tlist[0]) for tlist in e2t.values() for badt in tlist[1:] if len(tlist)>1)
    for p in pieces:
        filtered = set([tx_replace.get(t,t) for t in p.transcripts])
        transcript_ids |= filtered
        p.transcripts = list(filtered)
    transcript_ids = list(transcript_ids)
    gene_ids = list(set(e.gene_id for e in exons))


    #--- Remake the transcript-pieces mapping
    tp_map = {}
    for p in pieces:
        for tx in p.transcripts:
            tp_map.setdefault(tx,[]).append(p)
    #--- Remake the transcripts-exons mapping
    te_map = {}  # map {transcript: [exons]}
    for exon in exons:
        txs = exon.transcripts
        for t in txs:
            te_map.setdefault(t,[]).append(exon)


    #--- Get all reads from this chunk - iterator
    ckreads = sam.fetch(chrom, exons[0].start, lastend)


    #--- Count reads in each piece
    def count_reads(pieces,ckreads):
        #NH = [1.0/t[1] for t in read.tags if t[0]=='NH']+[1]
        try: read = ckreads.next()
        except StopIteration: return
        pos = read.pos
        rlen = read.rlen
        for p in pieces:
            if p.start + rlen < pos:
                continue
            while pos < p.end + rlen:
                p.count += 1.
                try: read = ckreads.next()
                except StopIteration: return
                pos = read.pos

    count_reads(pieces,ckreads)
    #--- Calculate RPK
    for p in pieces:
        p.rpk = toRPK(p.count,p.length)


    def estimate_expression(feat_class, pieces, ids):
        #--- Build the exons-transcripts structure matrix:
        # Lines are exons, columns are transcripts,
        # so that A[i,j]!=0 means "transcript Tj contains exon Ei".
        if feat_class == Gene:
            is_in = lambda p,g: g in p.gene_id.split('|')
        elif feat_class == Transcript:
            is_in = lambda p,t: t in p.transcripts
        n = len(pieces)
        m = len(ids)
        A = zeros((n,m))
        for i,p in enumerate(pieces):
            for j,f in enumerate(ids):
                A[i,j] = 1 if is_in(p,f) else 0
        #--- Build the exons scores vector
        E = asarray([p.rpk for p in pieces])
        #--- Solve for RPK
        T = lsqnonneg(A,E)
        #--- Store result in *feat_class* objects
        feats = []
        for i,f in enumerate(ids):
            exs = sorted([p for p in pieces if is_in(p,f)], key=lambda x:(x.start,x.end))
            flen = sum(p.length for p in pieces if is_in(p,f))
            feats.append(Transcript(name=f, start=exs[0].start, end=exs[-1].end,
                    length=flen, rpk=T[0][i], count=fromRPK(T[0][i],flen),
                    chrom=exs[0].chrom, gene_id=exs[0].gene_id, gene_name=exs[0].gene_name))
        return feats

    genes = estimate_expression(Gene, pieces, gene_ids)
    transcripts = estimate_expression(Transcript, pieces, transcript_ids)

    for t in transcripts: print t,t.count,t.rpk
    for g in genes: print g,g.count,g.rpk

    return genes,transcripts


def rnacount(bamname, annotname):
    """Annotation in GTF format, assumed to be sorted at least w.r.t. chrom name."""
    sam = pysam.Samfile(bamname, "rb")
    annot = open(annotname, "r")

    row = annot.readline().strip()
    exon0 = parse_gtf(row)
    chrom = exon0.chrom
    lastchrom = chrom

    while row:

        # Load all GTF exons of a chromosome in memory and sort
        chrexons = []
        while chrom == lastchrom:  # start <= lastend and
            exon = parse_gtf(row)
            if exon.end - exon.start > 1 :
                chrexons.append(exon)
            row = annot.readline().strip()
            if not row:
                break
        lastchrom = chrom
        chrexons.sort(key=lambda x: (x.start,x.end))
        print ">> Chromosome", chrom

        # Process chunks of overlapping exons / exons of the same gene
        lastend = chrexons[0].end
        lastgeneid = ''
        ckexons = []
        for exon in chrexons:
            # Store
            if (exon.start <= lastend) or (exon.gene_id == lastgeneid):
                ckexons.append(exon)
            # Process the stored chunk of exons
            else:
                process_chunk(ckexons, sam, chrom, lastend)
                ckexons = [exon]
            lastend = max(exon.end,lastend)
            lastgeneid = exon.gene_id
        process_chunk(ckexons, sam, chrom, lastend)

    annot.close()
    sam.close

######################################################################

bamname = "testfiles/gapdhKO.bam"
annotname = "testfiles/mm9_mini.gtf"

rnacount(bamname,annotname)



