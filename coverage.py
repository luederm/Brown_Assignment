import subprocess
import os
import luigi
from luigi import mock

class Read():
    '''
        Represents a sequencing read. Constructed from lines of fastq file.
    '''
    def __init__(self, seqID, seq, quality):
        self.seqID = seqID.trim()
        self.seq = seq.trim()
        self.quality = seq.trim()

    def phred33ToQ(self, char):
        ''' Convert a phred 33 encoded character to quality score '''
        return ord(char) - 33



class FilterReads(luigi.Task):
    '''
        Filter out reads with mean quality score less than Q
    '''
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()
    t = luigi.IntParameter()
    q = luigi.IntParameter()

    def requires(self):
        return []

    def output(self):
        return [luigi.LocalTarget('r1.filt.fastq'), luigi.LocalTarget('r2.filt.fastq')] #TODO: base filenames off r1 and r2 parameters


class Index(luigi.Task):
    '''
        Create an index for the reference genome (Using BWA and SAMtools)
    '''
    x = luigi.Parameter()

    def requires(self):
        return []

    def output(self):
        return [luigi.LocalTarget(str(self.x) + ft) for ft in ['.amb', '.ann', '.bwt', '.fai', '.pac', '.sa']]

    def run(self):
        print("Indexing reference with BWA and SAMtools...")
        subprocess.run(['bwa', 'index', '-a', 'bwtsw', self.x])
        subprocess.run(['samtools', 'faidx', self.x])


class AlignToRefGenome(luigi.Task):
    '''
        Use BWA-MEM to align reads to the reference genome
    '''
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()
    x = luigi.Parameter()
    t = luigi.IntParameter()

    def requires(self):
        return [Index(x = self.x)]

    def output(self):
        return luigi.LocalTarget("alignment.bam")

    def run(self):
        print("Aligning reads to reference genome with BWA-MEM...")
        bwaOut = subprocess.Popen(['bwa', 'mem', '-t', str(self.t), self.x, self.r1, self.r2], stdout=subprocess.PIPE)
        samToBam = subprocess.Popen(['samtools', 'view', '-b', '-o', 'alignment.bam'], stdin=bwaOut.stdout)
        samToBam.wait()


class SortBAM(luigi.Task):
    '''
        Sort BAM file using SAMtools
    '''
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()
    x = luigi.Parameter()
    t = luigi.IntParameter()

    def requires(self):
        return [AlignToRefGenome(r1 = self.r1, r2 = self.r2, x = self.x, t = self.t)]

    def output(self):
        return luigi.LocalTarget("alignment.sort.bam")

    def run(self):
        print("Sorting BAM file...")
        subprocess.run(['samtools', 'sort', '-o', 'alignment.sort.bam', '--threads', str(self.t), self.input()[0].fn])


class CalcCoverage(luigi.Task):
    '''
        Uses BAM file to calculate coverage for each bp in the reference genome.
        Generates coverage.tsv file: A tab-delimited file with columns “position” (in the reference sequence)
            and “count”. Count is defined as the number of mapped reads that contain the given
            reference position in their mapping.
        Compiles some stats/metrics for correlation file generation later in pipeline
    '''
    r1 = luigi.Parameter()
    r2 = luigi.Parameter()
    x = luigi.Parameter()
    t = luigi.IntParameter()

    def requires(self):
        return [SortBAM(r1 = self.r1, r2 = self.r2, x = self.x, t = self.t)]

    def output(self):
        return [luigi.LocalTarget("coverage.tsv"), luigi.LocalTarget("base_anova.csv"),
                luigi.LocalTarget("GC_vs_cov.csv")]

    def execute(self, cmd):
        popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
        for stdout_line in iter(popen.stdout.readline, ""):
            yield stdout_line
        popen.stdout.close()
        return_code = popen.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)

    def run(self):
        with self.output()[0].open('w') as coverage_tsv:
            coverage_tsv.truncate()
            perBaseCov = {'G': [], 'C': [], 'T': [], 'A': []}
            i = 1
            averageGC = []
            averageCov = []
            GC_sum = 0
            cov_sum = 0

            for baseDat in self.execute(["samtools", "mpileup", "-f", self.x, "-R", self.input()[0].fn]):
                baseDat = baseDat.split('\t')
                coverage_tsv.write("%s\t%s\n" % (baseDat[1], baseDat[3]))

                # Create a list of coverage values for each base type
                if str.upper(baseDat[2]) in ('G', 'C', 'T', 'A'):
                    perBaseCov[str.upper(baseDat[2])].append(int(baseDat[3]))

                # Count the number of G/C bases (for every 100 bases)
                if baseDat[2] in ('G', 'C'):
                    GC_sum += 1

                # Sum the coverage for each base (for every 100 bases)
                cov_sum += int(baseDat[3])

                # Save average GC% and coverage for every 100 bases
                if (i % 100 == 0):
                    averageGC.append(GC_sum / 100)
                    GC_sum = 0
                    averageCov.append(cov_sum / 100)
                    cov_sum = 0

                i += 1

            with self.output()[1].open('w') as base_stats:
                base_stats.write('"Base","Coverage"\n')
                for base in perBaseCov.keys():
                    for cov in perBaseCov[base]:
                        base_stats.write("%s,%i\n" % (base, cov))

            with self.output()[2].open('w') as gc_stats:
                gc_stats.write('"avg.GC%","avg.cov"\n')
                for i in range(0, len(averageCov)):
                    gc_stats.write("%f,%f\n" % (averageGC[i], averageCov[i]))

'''
# Example
os.chdir("/home/luederm/Desktop/Brown-Assignment")
coverage_tsv = open('coverage.tsv', 'w')
coverage_tsv.truncate()
totalCov = {'G':0, 'C':0, 'T':0, 'A':0}
numBases = {'G':0, 'C':0, 'T':0, 'A':0}
for baseDat in execute(["samtools", "mpileup", "-f", "Ref/ref.fasta", "-R", "alignment.sort.bam"]):
    baseDat = baseDat.split('\t')
    coverage_tsv.write("%s\t%s\n" % (baseDat[1], baseDat[3]))
    totalCov[baseDat[2]] += int(baseDat[3])
    numBases[baseDat[2]] += 1

for base in ['G', 'C', 'T', 'A']:
    print(base)
    print(totalCov[base])
    print(numBases[base])
    print(totalCov[base] / numBases[base])
    print('\n')
'''