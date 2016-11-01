# File: DNAFlexRBP.py 
#
# Copyright (C) 2016 Marco Pasi <marco.pasi@nottingham.ac.uk> 
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
#"DNAFlexRBP.py v0.1 (C) 2016 Marco Pasi"
#
#

from mug import datatypes as mug_datatypes
from DNAFlexibility import DNAFlexibility

from pycompss.api.task import task
from pycompss.api.constraint import constraint
from pycompss.api.parameter import IN, OUT

import urllib2, getopt
import numpy as np
from scipy.linalg import block_diag

#------------------------------------------------------------------------------
__doc__ = """\
Predict  the  flexibility  of  a  canonical  B-form  DNA  oligomer  of
arbitrary  sequence  by  using  the  NAFlex  webserver.   NAFlex  uses
MD-derived,  tetranucleotide-dependent  stiffness  matrices  (inverted
covariance  matrices) as  a  function of  the Curves+  inter-base-pair
helical parameters (i.e. Shift, Slide,  Rise, Tilt, Roll, Twist).  For
more  information  see the  following  references.  If you  find  this
software useful, please cite the following references.

Usage:
DNAFlexRBP.py [options] <SEQUENCE>

Options:
-h (--help)		Display this help.

References:
[1] Lavery,R., Moakher,M., Maddocks,J.H., Petkeviciute,D.  and
Zakrzewska,K. (2009) Conformational analysis of nucleic acids
revisited: Curves+.  Nucleic Acids Res., 37, 5917-5929.
[2] Hospital,A., Faustino,I., Collepardo-Guevara,R., Gonzalez,C.,
Gelpi,J.L. and Orozco,M. (2013) NAFlex: a web server for the study of
nucleic acid flexibility. Nucleic Acids Res., 41, W47-W55.
"""

#------------------------------------------------------------------------------
bases=list("ATCG")
_other = {
    "A":"T",
    "C":"G",
    "T":"A",
    "G":"C",
    "R":"Y",
    "Y":"R"
    }
def other(seq):
    return [_other[x] for x in seq]
def wcc(seq):
    return "".join(reversed(other(seq)))

#------------------------------------------------------------------------------
class DNAFlexRBP(DNAFlexibility):

    output_data_type = mug_datatypes.Matrix
    configuration = {}

    ABC_info = map(lambda x:x.split(","), """\
AAAA,AAAA,8;AAAC,AAAC,9;AAAG,GAAA,10;AAAT,AAAT,9;AACA,AAAC,10;AACC,TGGT,27;\
AACG,GAAC,10;AACT,TAAC,10;AAGA,GAAA,7;AAGC,CAAG,10;AAGG,AGGA,8;AAGT,TAAG,10;\
AATA,AAAT,10;AATC,CAAT,10;AATG,GAAT,10;AATT,TAAT,10;ACAA,AAAC,7;ACAC,TGTG,29;\
ACAG,TGTC,29;ACAT,TGTA,29;ACCA,TGGT,28;ACCC,GGGT,28;ACCG,CGGT,28;ACCT,AGGT,28;\
ACGA,GAAC,7;ACGC,CGTG,29;ACGG,CGGA,8;ACGT,CGTA,8;ACTA,TAAC,7;ACTC,AGTG,29;\
ACTG,AGTC,29;ACTT,TAAG,27;AGAA,GAAA,8;AGAC,TGTC,27;AGAG,AGAG,9;AGAT,GATA,8;\
AGCA,CAAG,7;AGCC,TGGC,27;AGCG,AGCG,9;AGCT,AGCT,9;AGGA,AGGA,9;AGGC,AGGC,9;\
AGGG,GGGA,8;AGGT,AGGT,9;AGTA,TAAG,7;AGTC,AGTC,9;AGTG,AGTG,9;AGTT,TAAC,27;\
ATAA,AAAT,7;ATAC,TGTA,27;ATAG,GATA,10;ATAT,TATA,8;ATCA,CAAT,7;ATCC,TGGA,27;\
ATCG,TCGA,8;ATCT,GATA,29;ATGA,GAAT,7;ATGC,ATGC,9;ATGG,TGGA,8;ATGT,TGTA,8;\
ATTA,TAAT,28;ATTC,GAAT,28;ATTG,CAAT,28;ATTT,AAAT,28;CAAA,AAAC,8;CAAC,TGGT,30;\
CAAG,CAAG,9;CAAT,CAAT,9;CACA,TGTG,28;CACC,GGGT,27;CACG,CGTG,28;CACT,AGTG,28;\
CAGA,TGTC,30;CAGC,TGGC,30;CAGG,AGGC,8;CAGT,AGTC,8;CATA,TGTA,30;CATC,TGGA,30;\
CATG,ATGC,8;CATT,GAAT,27;CCAA,TGGT,29;CCAC,GGGT,30;CCAG,TGGC,29;CCAT,TGGA,29;\
CCCA,GGGT,29;CCCC,GGGG,29;CCCG,GGGC,29;CCCT,GGGA,29;CCGA,CGGT,29;CCGC,GGGC,30;\
CCGG,CGGC,8;CCGT,CGGA,29;CCTA,AGGT,29;CCTC,GGGA,30;CCTG,AGGC,29;CCTT,AGGA,29;\
CGAA,GAAC,8;CGAC,CGGT,30;CGAG,AGCG,7;CGAT,TCGA,29;CGCA,CGTG,30;CGCC,GGGC,27;\
CGCG,CGCG,9;CGCT,AGCG,28;CGGA,CGGA,9;CGGC,CGGC,9;CGGG,GGGC,8;CGGT,CGGT,9;\
CGTA,CGTA,9;CGTC,CGGA,30;CGTG,CGTG,9;CGTT,GAAC,27;CTAA,TAAC,8;CTAC,AGGT,30;\
CTAG,AGCT,7;CTAT,GATA,27;CTCA,AGTG,30;CTCC,GGGA,27;CTCG,AGCG,30;CTCT,AGAG,28;\
CTGA,AGTC,30;CTGC,AGGC,30;CTGG,TGGC,8;CTGT,TGTC,8;CTTA,TAAG,28;CTTC,AGGA,30;\
CTTG,CAAG,28;CTTT,GAAA,27;GAAA,GAAA,9;GAAC,GAAC,9;GAAG,AGGA,7;GAAT,GAAT,9;\
GACA,TGTC,28;GACC,CGGT,27;GACG,CGGA,7;GACT,AGTC,28;GAGA,AGAG,8;GAGC,AGCG,8;\
GAGG,GGGA,7;GAGT,AGTG,8;GATA,GATA,9;GATC,TCGA,7;GATG,TGGA,7;GATT,CAAT,27;\
GCAA,CAAG,8;GCAC,CGTG,27;GCAG,AGGC,7;GCAT,ATGC,28;GCCA,TGGC,28;GCCC,GGGC,28;\
GCCG,CGGC,28;GCCT,AGGC,28;GCGA,AGCG,10;GCGC,CGCG,8;GCGG,GGGC,7;GCGT,CGTG,8;\
GCTA,AGCT,29;GCTC,AGCG,29;GCTG,TGGC,7;GCTT,CAAG,27;GGAA,AGGA,10;GGAC,CGGA,10;\
GGAG,GGGA,10;GGAT,TGGA,10;GGCA,AGGC,10;GGCC,CGGC,10;GGCG,GGGC,10;GGCT,TGGC,10;\
GGGA,GGGA,9;GGGC,GGGC,9;GGGG,GGGG,8;GGGT,GGGT,9;GGTA,AGGT,10;GGTC,CGGT,10;\
GGTG,GGGT,10;GGTT,TGGT,10;GTAA,TAAG,8;GTAC,CGTA,10;GTAG,AGGT,7;GTAT,TGTA,10;\
GTCA,AGTC,10;GTCC,CGGA,27;GTCG,CGGT,7;GTCT,TGTC,10;GTGA,AGTG,10;GTGC,CGTG,10;\
GTGG,GGGT,7;GTGT,TGTG,8;GTTA,TAAC,28;GTTC,GAAC,28;GTTG,TGGT,7;GTTT,AAAC,28;\
TAAA,AAAT,8;TAAC,TAAC,9;TAAG,TAAG,9;TAAT,TAAT,9;TACA,TGTA,28;TACC,AGGT,27;\
TACG,CGTA,28;TACT,TAAG,30;TAGA,GATA,7;TAGC,AGCT,8;TAGG,AGGT,8;TAGT,TAAC,30;\
TATA,TATA,9;TATC,GATA,28;TATG,TGTA,7;TATT,AAAT,27;TCAA,CAAT,8;TCAC,AGTG,27;\
TCAG,AGTC,7;TCAT,GAAT,30;TCCA,TGGA,28;TCCC,GGGA,28;TCCG,CGGA,28;TCCT,AGGA,28;\
TCGA,TCGA,9;TCGC,AGCG,27;TCGG,CGGT,8;TCGT,GAAC,30;TCTA,GATA,30;TCTC,AGAG,29;\
TCTG,TGTC,7;TCTT,GAAA,30;TGAA,GAAT,8;TGAC,AGTC,27;TGAG,AGTG,7;TGAT,CAAT,30;\
TGCA,ATGC,10;TGCC,AGGC,27;TGCG,CGTG,7;TGCT,CAAG,30;TGGA,TGGA,9;TGGC,TGGC,9;\
TGGG,GGGT,8;TGGT,TGGT,9;TGTA,TGTA,9;TGTC,TGTC,9;TGTG,TGTG,9;TGTT,AAAC,27;\
TTAA,TAAT,8;TTAC,TAAG,29;TTAG,TAAC,29;TTAT,AAAT,30;TTCA,GAAT,29;TTCC,AGGA,27;\
TTCG,GAAC,29;TTCT,GAAA,29;TTGA,CAAT,29;TTGC,CAAG,29;TTGG,TGGT,8;TTGT,AAAC,30;\
TTTA,AAAT,29;TTTC,GAAA,28;TTTG,AAAC,29;TTTT,AAAA,29""".split(";"))
    ABC_tetrads = dict(zip(zip(*ABC_info)[0],
                               [(x[1], int(x[2])) for x in ABC_info]))

    def _ABC_tetrad_position(self, tetrad):
        """
        Return the 4-letter name of the ABC oligomer where the specified
        tetrad is found, and its position within the oligomer's sequence.
        """
        return self.ABC_tetrads[tetrad]
    
    def _get_flexibility_http(self, tetrad):
        """
        Download the stiffness matrix relative to the specified tetrad from
        the NAFlex website, and return it as a numpy array.
        """
        oligomer, position = self._ABC_tetrad_position(tetrad)
        if position > 18:
            position = 36 - position
        fullseq = "GC"+oligomer[2:]+oligomer+oligomer+oligomer+"GC"
        dinucleotide = fullseq[position-1:position+1]
        dinucleotide = dinucleotide + wcc(dinucleotide)
        
        url = "http://mmb.irbbarcelona.org/NAFlex2/getFile.php?"+\
          "fileloc=../NAFlex2/NAFlex-Data/NAFlex_parmBSC1/"+\
          "NAFlex_mu{oligomer}/STIFFNESS/FORCE_CTES/{dinucleotide}.{position}.cte&type=curves".format(
            oligomer=oligomer,
            dinucleotide=dinucleotide.lower(),
            position=position)

        response = urllib2.urlopen(url)
        
        if response.getcode() is not 200 or \
          response.info().get("content-length") == '0':
            exc_str = "<{url}>"
            if response.geturl() != url:
                exc_str += " [redirected to <{url}>]".format(
                    url = response.geturl())
            exc_str += " replied {httpcode}\nHEADERS:\n{headers}" 
            headers = "\n".join("{}:{}".format(*item)
                                for item in response.info().dict.iteritems())
            raise Exception(exc_str.format(
                url = url,
                httpcode = response.getcode(),
                headers = headers))

        return np.loadtxt(response)

    def _get_flexibility_rest(self, tetrad):
        pass

    def _get_flexibility(self, tetrad):
        """Return the flexibility of the specified tetrad."""
        return self._get_flexibility_http(tetrad)

    @task(sequence = IN, stiffness_matrix = OUT)
    @constraint(AppSoftware="numpy,scipy")
    def run(self, sequence):
        """
        Build the full-oligomer stiffness matrix by positioning the 6x6 
        stiffness matrices for each tetranucleotide (obtained by calling 
        DNAFlexibility's _get_flexibilities()) on the diagonal.
        """
        stiffness_matrix = block_diag(*list(self._get_flexibilities(sequence)))
        return stiffness_matrix


#------------------------------------------------------------------------------
def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
    except getopt.error, msg:
        sys.stderr.write(msg+"\nSee --help\n")
        sys.exit(1)
        
    for o, a in opts:
        if o in ("-h", "--help"):
            print __doc__
            sys.exit(0)

    if len(args) < 1:
        print __doc__
        sys.exit(1)

    print DNAFlexRBP().run(args[0])

if __name__ == "__main__":
    main()    
    
