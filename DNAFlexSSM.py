# File: DNAFlexSSM.py 
#
# Copyright (C) 2016 Marco Pasi <marco.pasi@nottingham.ac.uk>, Charles A. Laughton <charles.laughton@nottingham.ac.uk>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
#"DNAFlexSSM.py v0.1 (C) 2016 Marco Pasi"
#
#

from mug import datatypes as mug_datatypes
from DNAFlexibility import DNAFlexibility

from pycompss.api.task import task
from pycompss.api.constraint import constraint
from pycompss.api.parameter import IN, OUT

#------------------------------------------------------------------------------
__doc__ = """\
Predict  the  flexibility  of  a  canonical  B-form  DNA  oligomer  of
arbitrary sequence  by using the  Stiffening Spring Model  (SSM).  The
SSM  uses  MD-derived  tetranucleotide variances  as  parameters,  but
accounts   for  the   conformational   frustration  between   adjacent
tetranucleotides,    which    makes   the    prediction    effectively
hexanucleotide-dependent.   For  more  information see  the  following
references.  If  you  find  this  software  useful,  please  cite  the
following references.

References:
[1] Shkurti,A., Goni,R., Andrio,P., Breitmoser,E., Bethune,I.,
Orozco,M. and Laughton,C.A. (2016) pyPcazip: A PCA-based toolkit for
compression and analysis of molecular simulation data. SoftwareX.
[2] Pasi,M., Triantafyliou,S., Shkurti,A. and Laughton,C.A. (2016) The
Sequence-Dependence of DNA Structure and Flexibility - Lessons From
The muABC Dataset. in preparation
"""

#------------------------------------------------------------------------------
def read_fasta(fname):
    """Get the sequence from a fasta-format file."""
    seq = ""
    with open(fname,'r') as f:
        f.readline() # skip header
        for line in f:
            seq = seq + line[:-1]
    return seq

def seq2yr(seq):
    """Convert a DNA sequence to the RY alphabet, where R=A,G; Y=C,T."""
    result = ""
    for i in range(len(seq)):
        if seq[i]=='A' or seq[i] == 'G':
            result = result + 'R'
        else:
            result = result + 'Y'
    return result


#------------------------------------------------------------------------------
class fdat(dict):
    """The file format of the SSM parameters."""
    def __init__(self, fname=None):
        if fname is not None:
            self.read(fname)

    def read(self, fname):
        ret = {}
        with open(fname,'r') as f:
            for line in f:
                words = line.split()
                self[words[0]] = float(words[1])

    def write(self, fname):
        with open(fname,'w') as f:
            for k,v in sdict.iteritems():
                f.write("{} {:f}\n".format(k, v))

#------------------------------------------------------------------------------
class DNAFlexSSM(DNAFlexibility):

    output_data_type = mug_datatypes.Vector
    """Configuration: the names of the two parameter files."""
    configuration = {
        'tet_score_file':  'tetvariance.dat',
        'comp_score_file': 'YRcompatibility.dat'
        }

    def __init__(self, configuration = {}):
        super(DNAFlexSSM, self).__init__(configuration)
        "Load parameter files"
        self.tet_scores = fdat(self.configuration['tet_score_file'])
        self.comp_scores= fdat(self.configuration['comp_score_file'])

    def _predict(self, hexad):
        """
        Predict the SSM variance of the specified hexad.
        """
        c1 = self.comp_scores[seq2yr(hexad[0:5])]
        c2 = self.comp_scores[seq2yr(hexad[1:6])]
        score = self.tet_scores[hexad[1:5]]/(1+(c1*c1+c2*c2)*0.5)
        return score
    
    def _get_flexibility(self, hexad):
        """Return the flexibility of the specified hexad."""
        return self._predict(hexad)

    @task(sequence = IN, flexibility = OUT)
    def run(self, sequence):
        """
        Build the full-oligomer flexibility vector by concatenating the 
        SSM variances of each tetranucleotide (obtained by calling 
        DNAFlexibility's _get_flexibilities()).
        """
        flexibility = list(self._get_flexibilities(sequence, k=6))
        return flexibility


#------------------------------------------------------------------------------
if __name__ == "__main__":
    pass
