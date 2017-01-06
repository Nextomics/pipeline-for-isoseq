 # -*- coding: UTF-8 -*-   
#!/usr/bin/env python
##huj@nextomics.org

import re
import os
import sys
import gzip
import signal
import logging
import argparse
import svgwrite
import subprocess
import time
import networkx as nx
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from itertools import combinations

START = 200
XBIN = 400
YBIN = 150
FLANK = 100
HEIGHT = 30
SPLICESITELENGTH = 6
CONVERTBIN = '/usr/bin/convert'
DONOR = ['GT', 'GC']
ACCEPTOR = ['AG']
#canonicalsplicesite = GT/AG, GC/AG, ATATC/AG, ATATC/AC, ATATC/AT, GTATC/AT or ATATC/AA
#intron au content

logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s] %(asctime)s %(message)s')

class TimedOutExc(Exception):
    pass


class HelpFormatter(argparse.RawDescriptionHelpFormatter,argparse.ArgumentDefaultsHelpFormatter):
    pass


class deco(object):
    @staticmethod
    def _deadline(timeout):
        def _decorate(f):
            def handler(signum, frame):
                 raise TimedOutExc
            def new_f(*args):
                signal.signal(signal.SIGALRM, handler)
                signal.alarm(timeout)
                try:
                    return f(*args)
                except TimedOutExc:
                    logging.warning('convent svg to png timeout, please convent the file\
                        [' +  args[0].outFile + '.svg] to png manually.')
                    return 0
                finally:
                    signal.alarm(0)
            return new_f
        return _decorate


class File(object):
    '''File Base Class'''

    def __init__(self, inFile):
        self.file = inFile
        self.baseName = os.path.basename(inFile)
        self.type = self._findType()
        self.sName = self._shortName()

    def _findType(self):
        self._fileExists()
        if self.baseName.endswith(('.gz', '.z', '.gzip', '.GZ', '.Z', '.GZIP')):
            return 'gz'
        else:
            return 'txt'

    def _shortName(self):
        lable = -1
        if self.type == 'gz':
            lable -= 1
        return '.'.join(self.baseName.split('.')[:lable]) + '.'

    def _openFile(self):
        if self.type == 'txt':
            return open(self.file)
        else:
            return gzip.GzipFile(self.file)

    def _closeFile(self, handle):
        handle.close()

    def _fileExists(self):
        if not os.path.isfile(self.file):
            logging.error('Error, ' + self.file + ' is not exists.')
            sys.exit(1)

    def readRecord(self):
        pass


class Table(File):
    '''table file'''
    def __init__(self, inFile):
        super(Table, self).__init__(inFile)
        self.records = self._readRecords()

    def _readRecord(self, field=0):
        handle = self._openFile()
        for line in handle:
            if line.startswith('#') or line == '\n':
                continue
            lines = line.strip().split('\t')
            yield lines[field]
        self._closeFile(handle)

    def _readRecords(self):
        records = []
        for name in self._readRecord():
            records.append(name)
        return set(records)


class Fasta(File):
    """fasta class"""
    def __init__(self, inFile):
        super(Fasta, self).__init__(inFile)
        self.records = self._readRecords()

    def _fastaLength(self):
        fastaNumber = fastaLength = 0
        for record in self.readRecord():
            fastaNumber += 1
            fastaLength += len(record[1])
        return (fastaNumber, fastaLength)

    def _cutFileNumber(self, number, outDir):
        outDirs = []
        lable = totalLength = 0
        fastaNumber, fastaLength = self._fastaLength()
        cutLength = fastaLength / number if fastaLength % number == 0 else fastaLength / number + 1

        for record in self.readRecord():
            if totalLength == 0 or totalLength >= cutLength:
                if totalLength != 0:
                    out.close()
                totalLength = 0
                lable += 1
                outFile = outDir + '/' + self.sName + str(lable) + '.cut'
                out = open(outFile, 'w')
                outDirs.append(outFile)
            print >>out, '>' + record[0] + '\n' + record[1]
            totalLength += len(record[1])
        else:
            out.close()
        return outDirs

    def _cutSeqNumber(self, number, outDir):
        outDirs = []
        lable = totalNumber = 0
        for record in self.readRecord():
            if totalNumber == 0 or totalNumber == number:
                if totalNumber != 0:
                    out.close()
                totalNumber = 0
                lable += 1
                outFile = outDir + '/' + self.sName + str(lable) + '.cut'
                out = open(outFile, 'w')
                outDirs.append(outFile)
            print >>out, '>' + record[0] + '\n' + record[1]
            totalNumber += 1
        else:
            out.close()
        return outDirs

    def lineFormat(self, seq, lineNumber=60):
        return '\n'.join([seq[i:i + lineNumber] for i in xrange(0, len(seq), lineNumber)])

    def _readRecord(self):
        fastaID = fastaSeq = ''
        handle = self._openFile()
        for line in handle:
            if line.startswith('>'):
                if fastaID:
                    yield (fastaID, str(Seq(fastaSeq).back_transcribe()).upper())
                fastaID, fastaSeq = line.strip().split()[0][1:], ''
            else:
                fastaSeq += line.strip()
        else:
            yield (fastaID, str(Seq(fastaSeq).back_transcribe()).upper())
        self._closeFile(handle)

    def _readRecords(self):
        records = {}
        for name, seq in self._readRecord():
            records[name] = seq
        return records

    def cutFile(self, number=20, mode=1, outDir=None):
        '''1: cut file with the specified number of subfile in total.
           2: cut file with the specified number of seqeunces in each sub file.
        '''
        if (not outDir):
            outDir = self.baseName + '.cut'
        if mode == 1:
            outDirs = self._cutFileNumber(number, outDir)
        elif mode == 2:
            outDirs = self._cutSeqNumber(number, outDir)
        else:
            loging.error('Error mode[' + str(mode) + '], only can be used with 1,2.')
            sys.exit(2)
        return outDirs


class Gtf(File):
    '''gtf class'''
    def __init__(self, inFile, regex=r'gene_id\s+([^;]+).*transcript_id\s+([^;]+)'):
        super(Gtf, self).__init__(inFile)
        self.regex = regex
        self.exonlength = dict()
        self.records = defaultdict(dict)
        self.geneIDs = defaultdict(dict)
        self.introns = defaultdict(dict)
        self.spliceSeq = defaultdict(dict)
        self.geneRange = defaultdict(list)

    def _reverse_complement(self, seq):
        return str(Seq(seq).reverse_complement())

    def getMaxLengthTranscript(self, cdnas):
        maxLength = 0
        transcript = ''
        for i in cdnas:
            if self.exonlength[i] > maxLength:
                maxLength = self.exonlength[i]
                transcript = i
        return transcript

    def parseGtf(self, feature='exon'):
        handle = self._openFile()
        for line in handle:
            geneID = transcriptID = ''
            lines = line.strip().split("\t")
            if line.startswith('#') or line == '\n' or lines[2].lower() != feature.lower():
                continue
            pattern = re.compile(self.regex, re.I)
            group = pattern.search(lines[8])
            if not group:
                logging.error('Error,escape the regex[' + self.regex + '], please check.')
                sys.exit(3)
            elif len(group.groups()) == 2:
                geneID, transcriptID = [i.strip('"') for i in group.group(1, 2)]
            elif len(group.groups()) == 1:
                geneID, transcriptID = '1', group.group(1)
            self.geneIDs[geneID][transcriptID] = lines[0]
            self.geneIDs[transcriptID][geneID] = lines[0]

            start, end = int(lines[3]), int(lines[4])
            start, end = (end, start) if start > end else (start, end)
            if transcriptID not in self.records[lines[0]]:
                self.exonlength[transcriptID] = 0
                self.records[lines[0]][transcriptID] = []
            self.exonlength[transcriptID] += end - start + 1
            self.records[lines[0]][transcriptID].append([start, end, lines[6]])
        self._closeFile(handle)

        # deleteSingleTranscrip = defaultdict(list) # we do not consider single exon transcript
        for scaffold in self.records:
            for transcriptID in self.records[scaffold]:
                self.records[scaffold][transcriptID] = sorted(self.records[scaffold][transcriptID], key=lambda x: x[0])

    def readIntron(self, args):
        '''only accept feature==exon/CDS'''
        for scaffold in self.records:
            for transcriptID in self.records[scaffold]:
                if len(self.records[scaffold][transcriptID]) != 1:
                    self.introns[scaffold][transcriptID] = []
                    sortCoor = self.records[scaffold][transcriptID]
                    i = 0
                    j = len(sortCoor) - 1
                    while (i < j): # cannot use for loop
                        if sortCoor[i + 1][0] - sortCoor[i][1] - 1 > args.minintronlength:
                            self.introns[scaffold][transcriptID].append([sortCoor[i][1] + 1, sortCoor[i + 1][0] - 1, sortCoor[0][2]])
                            i += 1
                        else:
                            logging.info('short intron & merge exon[' + str(sortCoor[i][0]) + ',' + str(sortCoor[i][1]) + '] with exon['\
                                 + str(sortCoor[i + 1][0]) + ',' +str(sortCoor[i + 1][1]) + '] in transcript[' + transcriptID + '].')
                            self.records[scaffold][transcriptID] = sortCoor[:i] + [[sortCoor[i][0], sortCoor[i+1][1],\
                                 sortCoor[i][2]]] + sortCoor[i+2:]
                            sortCoor = self.records[scaffold][transcriptID]
                            j -= 1
                            if len(sortCoor) == 1:
                                del self.introns[scaffold][transcriptID]
                                break

    def getSpliceSeq(self, genome):
        for scaffold in self.introns:
            for transcriptID in self.introns[scaffold]:
                coors = self.introns[scaffold][transcriptID]
                if coors:
                    strand = coors[0][2]
                    for start, end, strand in coors:
                        startSeq = genome[scaffold][start - 1 : start - 1 + SPLICESITELENGTH ]
                        endSeq = genome[scaffold][end - SPLICESITELENGTH : end]
                        if strand == "+":
                            self.spliceSeq[transcriptID][start] = [startSeq, 'donor']
                            self.spliceSeq[transcriptID][end] = [endSeq, 'acceptor']
                        else:
                            self.spliceSeq[transcriptID][start] = [self._reverse_complement(endSeq), 'acceptor']
                            self.spliceSeq[transcriptID][end] = [self._reverse_complement(startSeq), 'donor']

    def getGeneRange(self):
        outed = []
        for scaffold in self.records:
            for transcriptID in self.records[scaffold]:
                geneID = self.geneIDs[transcriptID].keys()[0]
                if geneID == '1':
                    self.geneRange[transcriptID] = [self.records[scaffold][transcriptID][0][0],
                                                    self.records[scaffold][transcriptID][-1][1],
                                                    self.records[scaffold][transcriptID][0][2]]
                elif geneID not in outed:
                    coors = Coordinate().getRange(sum([self.records[scaffold][transcriptID] for transcriptID in self.geneIDs[geneID]], []))
                    self.geneRange[geneID] = [coors[0], coors[-1], self.records[scaffold][transcriptID][0][2]]
                    outed.append(geneID)


class OutPut():
    def __init__(self, outDir):
        self.outDir = outDir + '/'
        self.outPicture = outDir + '/gene.cluster.picture' +  '/'
        self._checkDir(self.outPicture)

    def _checkDir(self, outDir):
        if not os.path.exists(outDir):
            os.makedirs(outDir)

    def outSpliceSeq(self, spliceSeq, geneIDs):
        acceptor = []
        donor = []
        temp = []
        for i in spliceSeq:
            chromosome = geneIDs[i].values()[0]
            for j in spliceSeq[i]:
                coor = chromosome + '-' + str(j)
                if coor in temp:
                    continue
                elif spliceSeq[i][j][1] == 'donor':
                    donor.append(spliceSeq[i][j][0])
                else:
                    acceptor.append(spliceSeq[i][j][0])
                temp.append(coor)
        return acceptor, donor

    def _outList(self, lists, sep='\n'):
        '''cannot include dict'''
        if lists and isinstance(lists[0], list):
            return sep.join([self._outList(i, "\t") for i in lists])
        else:
            return sep.join([str(i) for i in lists])

    def _outDict(self, dicts):
        '''can include list'''
        content = ''
        for i in dicts:
            if isinstance(dicts[i], dict):
                content += str(i) + '\t' + self._outDict(dicts[i])
            elif isinstance(dicts[i], list):
                content += str(i) + '\t' + self._outList(dicts[i], '\t')
            else:
                content += str(i) + '\t' + str(dicts[i])
            content += "\n"
        return content

    def _outDictDict(self, dicts):
        content = ''
        for i in dicts:
            for j in dicts[i]:
                content += str(i) + "\t" + str(j)  + "\t" + self._outList(dicts[i][j], '\t')
                content += "\n"
        return content

    def outAsCode(self, dicts, outFile):
        with open(self.outDir + outFile, 'w') as OUT:
            for i in dicts:
                for j in dicts[i]:
                    for con in dicts[i][j]:
                        print >>OUT, '\t'.join((con[0][1], 'Undefined', 'as_event',str(con[0][3]), str(con[0][4]), '.', con[0][2],\
                                 '.', 'transcript_id ' + i + ',' + j + '; gene_id ' + con[0][0] + '; structure ' + \
                                 con[1][0] + ',' + con[1][1] + '; splice_chain ' + con[1][2] + ',' + con[1][3]))

    def outAsType(self, cluster, outFile):
        with open(self.outDir + outFile, 'w') as OUT:
            for i in cluster.typeCodeNumber:
                print >>OUT, 'AS_Number\t' + i + '\t' + str(cluster.typeCodeNumber[i])
            print >>OUT, "\n###########################\n"

            for i in cluster.typeCodeGene:
                print >>OUT, 'Gene_Number\t' + i + '\t' + str(len(cluster.typeCodeGene[i]))
            print >>OUT, "\n###########################\n"

            for i in cluster.codeNumber:
                print >>OUT, 'Code_Number\t' + i + '\t' + str(cluster.codeNumber[i])

    def outFile(self, content, outFile):
        if outFile.upper() == "STDOUT":
            print content
        else: 
            with open(self.outDir + outFile, 'w') as OUT:
                print >>OUT, content


class SVG():
    def __init__(self, outFile):
        self.start = START
        self.xbin = XBIN
        self.flank = FLANK
        self.height = YBIN
        self.convert = CONVERTBIN
        self.max = 0
        self.min = 0
        self.factor = 1
        self.windows = 10
        self.outFile = outFile
        self.dwg = svgwrite.Drawing(filename=self.outFile + '.svg')
        self.markerPath = self.dwg.path(d="M0 0 L7 2 L6 4 L12 2 L6 0 L7 2", fill='black')
        self.leftMarker = self.setMarker()
        self.rightMarker = self.setMarker('0')

    def setMaxMin(self, coors):
        self.min = coors[0]
        self.max = coors[-1]

    def setFactor(self, factor=1):
        self.factor = (self.max - self.min)/5000.0 if self.max - self.min > 5000 and factor == 1 else factor

    def setMarker(self, orient="-180"):
        marker = self.dwg.marker(insert=(7, 2), size=(15, 15), orient=orient)
        marker.add(self.markerPath)
        self.dwg.defs.add(marker)
        return marker

    def setTranscript(self, transcriptID, coors, rectColour='#669933'):
        strand = coors[0][2]
        totalLength = coors[-1][-2] - coors[0][0] + 1
        line = self.dwg.add(self.dwg.line(
                                    start=(self.start - self.flank, self.height),
                                    end=((totalLength + coors[0][0] - self.min - 1)/self.factor +  self.start + self.flank, self.height),
                                    style="stroke-width:10", 
                                    fill="black",
                                    stroke='black'
                                    ))
        if strand == "+":
            line['marker-end'] = self.rightMarker.get_funciri()
        else:
            line['marker-start'] = self.leftMarker.get_funciri()

        self.dwg.add(self.dwg.text(
                                transcriptID, 
                                insert=((totalLength + coors[0][0] - self.min - 1)/self.factor + self.start + self.flank + 50 , self.height + 15), 
                                style="font-size:60px; font-family:Arial"
                                ))
        for s, e, st in coors:
            self.dwg.add(self.dwg.rect(
                        insert=((s - self.min)/self.factor + self.start, self.height - HEIGHT),
                        size=((e - s + 1)/self.factor, HEIGHT * 2 ),
                        ry=5,
                        fill=rectColour,
                        stroke='#999999',
                        stroke_width=3
                        ))
        self.height += YBIN
        
    def setGraduation(self):
        self.height += YBIN
        totalLength = (self.max - self.min + 1)/self.factor

        self.xbin = int(totalLength/self.windows)
        self.max = ((int(totalLength) - 1)/self.xbin + 1) * self.xbin + self.start + 1
        # self.max = int(self.start + totalLength + 1)

        self.dwg.add(self.dwg.line(
                        start=(self.start, self.height),
                        end=(self.max, self.height),
                        style = "stroke-width:10",
                        fill="black",
                        stroke='black'
                        ))
        for i in xrange(self.start, self.max, self.xbin):
            self.dwg.add(self.dwg.line(
                                start=(i, self.height),
                                end=(i, self.height - 20),
                                style = "stroke-width:10",
                                fill="black",
                                stroke='black'
                                ))
            if self.factor != 1 or (i - self.start)/self.xbin%2 == 0:
                self.dwg.add(self.dwg.text(
                                    int((i - self.start) * self.factor) + self.min,
                                    insert=(i - 60 , self.height + 60),
                                    style = "font-size:60px; font-family:Arial"
                                    ))
        
    def setGene(self, geneID):
        self.height += YBIN * 2
        self.dwg.add(self.dwg.text(
                            "Gene ID:" + geneID,
                            insert=(self.max/2, self.height),
                            style = "font-size:100px; font-family:Arial"
                            ))

    def save(self):
        self.dwg.save()
    
    # @deco._deadline(200)
    def setPng(self, resize="30%"):
        proc = subprocess.Popen(
            self.convert + ' -resize ' + resize + ' ' + self.outFile + '.svg' + ' ' + self.outFile + '.png',
            shell=True
            )
        poll_seconds = 2
        deadline = time.time() + 200
        while time.time() < deadline and proc.poll() == None:
            time.sleep(poll_seconds)
        if proc.poll() == None:
            proc.kill()
            logging.warning('convent svg to png timeout, please convent the file [' +  self.outFile + '.svg] to png manually.')
        stdoutdata, stderrdata = proc.communicate()
        # os.popen(self.convert + ' -resize ' + resize + ' ' + self.outFile + '.svg' + ' ' + self.outFile + '.png')


class Coordinate():
    def coorCmp(self, s1, s2, lable=0):
        '''
        Accept:     start1,end1,start2,end2,lable
        Lable = 0:  return overlapLength
        Lable = 1:  return overlapLength,overlapType
        lable = 2:  return overlapCoor,overlapType

        '''
        overlapLength = 0
        overlapType = 0
        overlapCoor = []
        (s1, e1, s2, e2) = (float(i) for i in (s1 + s2))
        assert s1 <= e1 and s2 <= e2, 'Error: need s1[%f] <= e1[%f] and s2[%f] <= e2[%f]' % (s1, e1, s2, e2)

        if e2 < s1:
            overlapType = 1
        elif s2 <= s1 and  s1 <= e2 <= e1:
            overlapLength = e2 - s1 + 1
            overlapType = 2
            overlapCoor = [s1, e2]
        elif s2 <= s1 and e2 >= e1:
            overlapLength = e1 - s1 + 1
            overlapType = 3
            overlapCoor = [s1, e1]
        elif s2 >= s1 and e2 <= e1:
            overlapLength = e2 - s2 + 1
            overlapType = 4
            overlapCoor = [s2, e2]
        elif  s1 <= s2 <= e1 and e2 >= e1:
            overlapLength = e1 - s2 + 1
            overlapType = 5
            overlapCoor = [s2, e1]
        elif s2 > e1:
            overlapType = 6
        else:
            logging.error('Error,connot find overlap type %f,%f,%f,%f' % (s1, e1, s2, e2))
            sys.exit(1)

        if lable == 0:
            return overlapLength
        elif lable == 1:
            return (overlapLength, overlapType)
        elif lable == 2:
            return (overlapCoor, overlapType)
        else:
            logging.error('Error,cannot accept the lable [' + str(lable) + '].')
            sys.exit(2)

    def getRange(self, coors):
        '''[[1,2],[3,4]]'''
        if isinstance(coors[0], list):
            sortCoors = sorted(coors, key=lambda x: x[0])
            return (sortCoors[0][0], sortCoors[-1][1])
        else:
            sortCoors = sorted(coors)
            return (sortCoors[0], sortCoors[-1])


class Cluster():
    def __init__(self):
        self.novelGene = list()
        self.noverCdna = list()
        self.asCode = defaultdict(dict)
        self.cluster = defaultdict(list)
        self.errorOrient = defaultdict(list)
        self.alternativeSplice = defaultdict(list)
        self.provedTranscript = defaultdict(dict)
        self.unprovedTranscript = defaultdict(list)
        self.typeCodeGene = defaultdict(list)
        self.codeNumber = dict()
        self.typeCodeNumber = {'ExonS':0, 'IntronR':0, 'AltD':0, 'AltA':0, 'AltP':0, 'Other':0}

    def setErrorOrient(self, gene1, gene2):
        self.errorOrient[gene1].append(gene2)

    def setCluster(self, gene, transcript):
        if isinstance(transcript, list):
            self.cluster[gene] += transcript
        else:
            self.cluster[gene].append(transcript)

    def setNovelGene(self, clusterGenes, genes):
        for gene in genes:
            if gene not in clusterGenes:
                self.novelGene.append(gene)

    def setNovelCdna(self, clusterCDNAs, cDNAs):
        for cDNA in cDNAs:
            if cDNA not in clusterCDNAs:
                self.noverCdna.append(cDNA)

    def _checkErrorOrient(self, gene1, gene2):
        return True if [geneID, gene1] in self.errorOrient or [geneID, gene2] in self.errorOrient else False

    def _extendOverlap(self, cooris, coorjs, rangeis, rangejs):
        i = -1
        for i in xrange(min(len(cooris), len(coorjs))):
            s1, e1, t1 = cooris[i]
            s2, e2, t2 = coorjs[i]
            overlapLength = Coordinate().coorCmp((s1, e1), (s2, e2))
            if overlapLength/(e1 - s1 + 1) < args.coverage or overlapLength/(e2 - s2 + 1) < args.coverage:
                return False
        i += 1
        if i == len(cooris) == len(coorjs):
            return True
        elif i == len(cooris):
            s, e, t = coorjs[i]
            if Coordinate().coorCmp((s, e), rangeis, 1)[1] in [3, 4]:
                return False
            else:
                return True
        else:
            s, e, t = cooris[i]
            if Coordinate().coorCmp((s, e), rangejs, 1)[1] in [3, 4]:
                return False
            else:
                return True

    def _cluster(self, coors):
        temp = []
        ps = pe = 0
        clusters = defaultdict(list)
        coors = sorted(coors, key=lambda x: x[0])
        for s, e, t in coors:
            if ps == 0:
                ps, pe = s, e
            elif ps <= s  < e <= pe:
                pass
            elif ps <= s <= pe <= e:
                pe = e
            elif pe < s:
                clusters[str([ps, pe])] = temp
                ps, pe = s, e
                temp = []
            else:
                logging.error('Error,connot find overlap type %f,%f,%f,%f' % (ps, pe, s, e))
                continue
            temp.append([s, e, t])
        clusters[str([ps, pe])] = temp
        return clusters

    def _sameTranscript(self, cooris, coorjs, args, rangeis, rangejs):
        for i in xrange(len(cooris)):
            s1, e1, t1 = cooris[i]
            for j in xrange(len(coorjs)):
                s2, e2, t2 = coorjs[j]
                overlapLength = Coordinate().coorCmp((s1, e1), (s2, e2))
                if overlapLength/(e1 - s1 + 1) > args.coverage and overlapLength/(e2 - s2 + 1) > args.coverage:
                    return self._extendOverlap(cooris[i + 1:], coorjs[j + 1:], rangeis, rangejs)
                elif overlapLength/(e1 - s1 + 1) > 1 - args.coverage or overlapLength/(e2 - s2 + 1) > 1 - args.coverage:
                    return False
                elif e1 < s2 and Coordinate().coorCmp((s1, e1), rangejs, 1)[1] in [3, 4]:
                    return False
                elif s1 > e2 and Coordinate().coorCmp((s2, e2), rangeis, 1)[1] in [3, 4]:
                    return False
        return False

    def _setprovedTranscript(self, genes, cdnas, geneID):
        if genes and cdnas:
            for i in genes:
                self.provedTranscript[geneID][i] = cdnas

    def _setunprovedTranscript(self, geneID, gene):
        for i in self.cluster[geneID]:
            if i not in self.provedTranscript[geneID] and i in gene:
                self.unprovedTranscript[geneID].append(i)

    def _setSplice(self, geneID, cDNA, deletecDNAs):
        '''get the max length transcript to represent the alternative splice transcripts'''
        deleteTranscripts = []
        for cdnas, genes in deletecDNAs:
            maxLengthTranscript = ''
            if not genes:
                maxLengthTranscript = cDNA.getMaxLengthTranscript(cdnas)
            [deleteTranscripts.append(i) for i in cdnas if i != maxLengthTranscript]
        for i in self.cluster[geneID]:
            if i not in deleteTranscripts and (geneID not in self.errorOrient or i not in self.errorOrient[geneID]):
                self.alternativeSplice[geneID].append(i)

    def _getCoor(self, i, gene, cDNA, chromosome):
        if i in gene.introns[chromosome]:
             return (gene.introns[chromosome][i], gene.spliceSeq[i], gene.records[chromosome][i])
        elif i in cDNA.introns[chromosome]:
            return (cDNA.introns[chromosome][i], cDNA.spliceSeq[i], cDNA.records[chromosome][i])
        elif i in gene.records[chromosome]:
            return (gene.records[chromosome][i], '', gene.records[chromosome][i])
        else:
            return (cDNA.records[chromosome][i], '', cDNA.records[chromosome][i])

    def setSplice(self, gene, cDNA, args):
        for geneID in self.cluster:
            G = nx.Graph()
            chromosome = gene.geneIDs[geneID].values()[0]
            for i, j in combinations(self.cluster[geneID], 2):
                coorintronis, spliceSeqis, coorexonis = self._getCoor(i, gene, cDNA, chromosome)
                coorintronjs, spliceSeqjs, coorexonjs = self._getCoor(j, gene, cDNA, chromosome)
                rangeis, rangejs = (coorexonis[0][0], coorexonis[-1][1]), (coorexonjs[0][0], coorexonjs[-1][1])
                if coorintronis[0][2] == coorintronjs[0][2] and self._sameTranscript(coorintronis,\
                        coorintronjs, args, rangeis, rangejs) and ((spliceSeqis and spliceSeqjs) or not (spliceSeqis and spliceSeqjs)):
                    G.add_edge(i, j)
            deletecDNAs = []
            familyTranscript = [g.nodes() for g in nx.connected_component_subgraphs(G)]
            for transcripts in familyTranscript:
                genes = []
                cdnas = []
                for transcript in transcripts:
                    if transcript in gene.records[chromosome]:
                        genes.append(transcript)
                    else:
                        cdnas.append(transcript)
                deletecDNAs.append([cdnas, genes])
                self._setprovedTranscript(genes, cdnas, geneID)
            self._setunprovedTranscript(geneID, gene.geneIDs[geneID])
            self._setSplice(geneID, cDNA, deletecDNAs)###add genes
            
    def setAsCode(self, gene, cDNA, args):
        for geneID in self.alternativeSplice:
            chromosome = gene.geneIDs[geneID].values()[0]
            for i ,j in combinations(self.alternativeSplice[geneID], 2):
                coorintronis, spliceSeqis, coorexonis = self._getCoor(i, gene, cDNA, chromosome)
                coorintronjs, spliceSeqjs, coorexonjs = self._getCoor(j, gene, cDNA, chromosome)
                rangeis, rangejs = (coorexonis[0][0], coorexonis[-1][1]), (coorexonjs[0][0], coorexonjs[-1][1])

                if coorintronis[0][2] != coorintronjs[0][2]:
                    logging.error(i  + ' & ' + j + ' orientation conflict.' )
                    continue

                elif spliceSeqis and spliceSeqjs:
                    strand = coorintronis[0][2]
                    cooris = sorted([k for k in sum(coorintronis, []) if k not in ['+', '-']])
                    coorjs = sorted([k for k in sum(coorintronjs, []) if k not in ['+', '-']])
                    coorintrons = coorintronis + coorintronjs
                    cluster = self._cluster(coorintrons)
                    if strand == "-":
                        cooris.reverse()
                        coorjs.reverse()
                    for k in cluster:
                        s, e = eval(k)
                        if Coordinate().coorCmp((s, e), rangeis, 1)[1] in [3, 4] and Coordinate().coorCmp((s, e), rangejs, 1)[1] in [3, 4]:
                            coors = sum(cluster[k], [])
                            coors = sorted([k for k in coors if coors.count(k) == 1 and k not in ['+', '-']])
                            if coors:
                                if strand == "-":
                                    coors.reverse()
                                if j not in self.asCode[i]:
                                    self.asCode[i][j] = []
                                self.asCode[i][j].append(((geneID, chromosome, strand, s, e), \
                                    self._getCode(coors, strand, cooris, spliceSeqis, coorjs, spliceSeqjs, args)))

    def _getCode(self, coors, strand, cooris, spliceSeqis, coorjs, spliceSeqjs, args):
        lable = codei = codej = codeCooi = codeCooj = ''
        for i in coors:
            indes = coors.index(i) + 1
            if i in cooris:
                lable = '^' if cooris.index(i)%2 == 0 else '-'
                codei += str(indes) + lable
                codeCooi += str(i) + lable
            elif i in coorjs:
                lable = '^' if coorjs.index(i)%2 == 0 else '-'
                codej += str(indes) + lable
                codeCooj += str(i) + lable
            else:
                logging.error('Position ' + str(i) + ' cannot find in ' + str(cooris) + ';' + str(coorjs) + '!')
        if not codei:
            codei = codeCooi = '0'
        if not codej:
            codej = codeCooj = '0'
        return (codei, codej, codeCooi, codeCooj)

    def typeAsCode(self, gene, cDNA, args):
        temp = []
        baseTranscript = []
        if os.path.exists(args.astypestat):
            baseTranscript = Table(args.astypestat).records
        for i in self.asCode:
            for j in self.asCode[i]:
                if not baseTranscript or i in baseTranscript or j in baseTranscript:
                    for code in self.asCode[i][j]:
                        geneID, chromosome, strand = code[0][:3]
                        codei, codej, codeCooi, codeCooj = code[1]
                        if chromosome + '#' + codeCooi + '#' + codeCooj in temp:
                            continue

                        if codei == '0':
                            if codej.endswith('-'):
                                self.typeCodeNumber['IntronR'] += 1
                                if geneID not in self.typeCodeGene['IntronR']:
                                    self.typeCodeGene['IntronR'].append(geneID)
                            else:
                                self.typeCodeNumber['ExonS'] += 1
                                if geneID not in self.typeCodeGene['ExonS']:
                                    self.typeCodeGene['ExonS'].append(geneID)
                        elif codej == '0':
                            if codei.endswith('-'):
                                self.typeCodeNumber['IntronR'] += 1
                                if geneID not in self.typeCodeGene['IntronR']:
                                    self.typeCodeGene['IntronR'].append(geneID)
                            else:
                                self.typeCodeNumber['ExonS'] += 1
                                if geneID not in self.typeCodeGene['ExonS']:
                                    self.typeCodeGene['ExonS'].append(geneID)
                        elif (codei == '1^' and codej == '2^') or (codei == '2^' and codej == '1^'):
                            self.typeCodeNumber['AltD'] += 1
                            if geneID not in self.typeCodeGene['AltD']:
                                self.typeCodeGene['AltD'].append(geneID)
                        elif (codei == '1-' and codej == '2-') or (codei == '2-' and codej == '1-'):
                            self.typeCodeNumber['AltA'] += 1
                            if geneID not in self.typeCodeGene['AltA']:
                                self.typeCodeGene['AltA'].append(geneID)
                        elif re.search(r'^([1234][\^-][1234][\^-])$', codei) and re.search(r'^([1234][\^-][1234][\^-])$', codej):
                            self.typeCodeNumber['AltP'] += 1
                            if geneID not in self.typeCodeGene['AltP']:
                                self.typeCodeGene['AltP'].append(geneID)
                        else:
                            self.typeCodeNumber['Other'] += 1
                            if geneID not in self.typeCodeGene['Other']:
                                self.typeCodeGene['Other'].append(geneID)

                        temp += [chromosome + '#' + codeCooi + '#' + codeCooj, chromosome + '#' + codeCooj + '#' + codeCooi]

                        if codei + ',' + codej not in self.codeNumber and codej + ',' + codei not in self.codeNumber:
                            self.codeNumber[codei + ',' + codej] = 0
                        if codei + ',' + codej in self.codeNumber:
                            self.codeNumber[codei + ',' + codej] += 1
                        else:
                            self.codeNumber[codej + ',' + codei] += 1


def main(args):

    coor = Coordinate()
    cluster = Cluster()
    out = OutPut(args.outdir)
    logging.info('Reading the genome fasta file...')
    fasta = Fasta(args.fasta)

    logging.info('Reading the gene file...')
    gene = Gtf(args.gtf)
    gene.parseGtf()
    gene.readIntron(args)
    gene.getGeneRange()
    gene.getSpliceSeq(fasta.records)

    logging.info('Reading the cDNA_match file...')
    cDNA = Gtf(args.input, r'ID=([^;]+)')
    cDNA.parseGtf('cDNA_match')
    cDNA.readIntron(args)
    cDNA.getGeneRange()
    cDNA.getSpliceSeq(fasta.records)

    for geneID in gene.geneRange:
        lable = True
        chromosome = gene.geneIDs[geneID].values()[0]
        for transcriptID in cDNA.records[chromosome]:
            if coor.coorCmp(gene.geneRange[geneID][:2], cDNA.geneRange[transcriptID][:2]):
                cluster.setCluster(geneID, transcriptID)
                if lable:
                    cluster.setCluster(geneID, gene.geneIDs[geneID].keys())
                    lable = False
                if gene.geneRange[geneID][-1]  != cDNA.geneRange[transcriptID][-1]:
                    cluster.setErrorOrient(geneID, transcriptID)
    cluster.setNovelGene(cluster.cluster.keys(), gene.geneRange.keys())
    cluster.setNovelCdna(sum(cluster.cluster.values(), []), cDNA.geneRange.keys())
    cluster.setSplice(gene, cDNA, args)

    out.outFile(out._outDict(cluster.cluster), 'transcript.cluster.list.txt')
    out.outFile(out._outList(cluster.novelGene), 'unproved.gene.list.txt')
    out.outFile(out._outList(cluster.noverCdna), 'novel.gene.list.txt')
    out.outFile(out._outDict(cluster.errorOrient), 'error.orient.list.txt')
    out.outFile(out._outDict(cluster.alternativeSplice), 'alternative.splice.list.txt')
    out.outFile(out._outDictDict(cluster.provedTranscript), 'proved.transcript.list.txt')
    out.outFile(out._outDict(cluster.unprovedTranscript), 'unproved.transcript.list.txt')

    if args.outsplice:
        acceptor, donor = out.outSpliceSeq(dict(gene.spliceSeq, **cDNA.spliceSeq), dict(gene.geneIDs, **cDNA.geneIDs))
        out.outFile(out._outList(acceptor), 'acceptor.list.txt')
        out.outFile(out._outList(donor), 'donor.list.txt')

    if args.outpicture:
        for geneID in cluster.cluster:
            svg = SVG(out.outPicture + geneID)
            svg.setMaxMin(coor.getRange(sum([cDNA.geneRange[transcriptID][:2] for transcriptID in \
                            cluster.cluster[geneID] if transcriptID not in gene.geneIDs], gene.geneRange[geneID][:2])))
            svg.setFactor()
            for transcriptID in gene.geneIDs[geneID]:
                chromosome = gene.geneIDs[geneID][transcriptID]
                coors = gene.records[chromosome][transcriptID]
                svg.setTranscript(transcriptID, coors)
            for transcriptID in cluster.cluster[geneID]:
                if transcriptID not in gene.geneIDs[geneID]:
                    chromosome = cDNA.geneIDs[transcriptID]['1']
                    coors = cDNA.records[chromosome][transcriptID]
                    svg.setTranscript(transcriptID, coors, 'red')
            svg.setGraduation()
            svg.setGene(geneID)
            svg.save()
            svg.setPng()

    if args.ascode:
        cluster.setAsCode(gene, cDNA, args)
        out.outAsCode(cluster.asCode, 'splice.ascode.list.txt')
    if args.astypestat:
        cluster.typeAsCode(gene, cDNA, args)
        out.outAsType(cluster, 'splice.ascode.stat.txt')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class = HelpFormatter,
        description = '''
name:
    alternative_splice.py -- Alternative Splicing Analysis Package


description:
    Alternative Splice Code:

        Each AS event is assigned an AS code according to the 
        relative position of the alternative splice sites that 
        are involved in the splicing variation.A variation in the
        exon-intron structure is detected by comparing all overlapping
        transcripts in a pairwise fashion. For a given pair, a splice
        site is said alternative if it is not used in both transcripts.
        In the same transcript comparison, considering all splice sites
        sorted by their genomic position, an AS event is defined as a
        maximal succession of alternative splice sites. Then, each 
        alternative splice site is assigned a number (1, 2, 3, ..) 
        according to its relative position in the event, and a symbol
        that depends on its type. To denote a donor site a "^" sign 
        is used (depicting the spliced out intron downstream of the
        donor site) and a "-" symbol characterizes an acceptor site.
        A splice site denoted "3-" therefore designates the 3rd splice
        site in an event, which is an acceptor. The AS code is built
        by writting down the splice sites of each transcript, separated
        by a comma. A 0 is used to denote transcripts that do not involve
        any of the alternative splice sites (e.g., in case of an exon 
        skipping event).

    Alternative Splice Type Statistics:

        The coordinates of reliable introns and exons were com-
        pared in a pairwise fashion in order to identify candidates
        for AS events. For intron/intron comparison, if two
        introns had the same 3'-end but a different 5'-end, this
        event was classified as AltD. If two introns differed only in
        the 3'-ends, this event was classified as AltA. AltP events
        refer to introns overlapping with each other but with both
        5'- and 3'-ends differing. For intron/exon comparisons, if
        an intron was completely covered by an exon, the event
        was classified as IntronR. If an exon was completely cov-
        ered by an intron, the event was classified as ExonS.


references
    [1] Sammeth M, Foissac S, Guigó R. A general definition and 
        nomenclature for alternative splicing events[J]. PLoS 
        Comput Biol, 2008, 4(8): e1000147.
    [2] Wang B B, Brendel V. Genomewide comparative analysis of
        alternative splicing in plants[J]. Proceedings of the 
        National Academy of Sciences, 2006, 103(18): 7175-7180.


demo
    alternative_splice.py -i cDNA.align -g gene.gtf -f genome.fa -os -op -as -ats
    '''
    )
    parser.add_argument('-i','--input',metavar = 'FILE', required = True,
        help = 'set the align file with GFF3 cDNA_match format.')
    parser.add_argument('-g','--gtf',metavar = 'FILE', required = True,
        help = 'set the gene file with gtf format.')
    parser.add_argument('-f','--fasta',metavar = 'FILE', required = True,
        help = 'set the genome file with fasta format.')

    parser.add_argument('-as','--ascode',action='store_true',
        help = 'caculate the alternative splice code.')
    parser.add_argument('-ats','--astypestat',metavar = 'FILE|BOOL',
        help = 'output alternative splice type statistics with given constitutive splicing ID table or For all.')

    parser.add_argument('-op','--outpicture',action='store_true',
        help = 'output the alternative splice picture svg & png.')
    parser.add_argument('-os','--outsplice',action='store_true',
        help = 'output all splice site sequences.')
    parser.add_argument('-ca','--canonical',action='store_true',
        help = 'consider only introns with canonical splice sites.')
    parser.add_argument('-c','--coverage',metavar = 'float',type = float,default = 0.9,
        help = 'filter exon if overlap with others and coverage >= $FILTER.')#, only for 5/3-end or non canonical splice sites.')

    parser.add_argument('-mil', '--minintronlength', type=int, default = 9,
        help = 'set the minimum length for one internal intron.')
    parser.add_argument('-t','--type',metavar = 'exon',default = 'exon',choices = ['exon', 'CDS'],
        help = 'set the accepted type.')
    parser.add_argument('-o','--outdir',metavar = 'DIRECTORY',default = './alternative.splice.result',
        help = 'set the output directory.')
    args = parser.parse_args()
    main(args)
