# File: DNAFlexRB.py 
#
# Copyright (C) 2016 Marco Pasi <marco.pasi@nottingham.ac.uk> 
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
#"DNAFlexRB.py v0.1 (C) 2016 Marco Pasi"
#
#

from mug import datatypes as mug_datatypes
from DNAFlexibility import DNAFlexibility

from pycompss.api.task import task
from pycompss.api.constraint import constraint
from pycompss.api.parameter import IN, OUT

from subprocess import Popen, PIPE
import tempfile, sys
import numpy as np

#------------------------------------------------------------------------------
__doc__ = """\
Predict  the  flexibility  of  a  canonical  B-form  DNA  oligomer  of
arbitrary sequence by using the  cgDNA tool.  cgDNA uses an MD-derived
parameter  set to  predict stiffness  matrices  as a  function of  the
Curves+ intra-base-pair  and inter-base-pair helical  parameters.  For
more  information see  the  following references.   If  you find  this
software useful, please cite the following references. Requires a 
working installation of GNU Octave, and of cgDNA in the Octave path.

Usage:
DNAFlexRB.py [options] <SEQUENCE>

Options:
-h (--help)		Display this help.

References:
[1] Gonzalez,O., Petkeviciute,D.  and Maddocks,J.H.  (2013) A
sequence-dependent rigid-base model of DNA. J Chem Phys, 138, 055102.
[2] Petkeviciute,D., Pasi,M., Gonzalez,O., and Maddocks,J.H. (2014)
cgDNA: a software package for the prediction of sequence-dependent
coarse-grain free energies of B-form DNA. Nucleic Acids Res., 42, e153-
"""

#------------------------------------------------------------------------------
class DNAFlexRB(DNAFlexibility):

    output_data_type = mug_datatypes.Matrix
    """Configuration: the name of the parameter file."""
    configuration = {
        'cgDNApath': None,
        'cgDNAparamset': 'cgDNAparamset2.mat'
        }
    
    def _is_error(self, stderr):
        """
        Parse the standard error stream from an Octave session to assess 
        whether an error was reported.
        """
        for line in stderr.split("\n"):
            if line.startswith("error:"): return True
        return False
        
    def _get_stiffness_matrix(self, sequence):
        """
        Generate the full-oligomer stiffness matrix by running an octave
        script that uses the cgDNA code.
        """

        with tempfile.NamedTemporaryFile(mode="w", suffix=".m") as octave_script,\
            tempfile.NamedTemporaryFile(mode="w", suffix=".dat") as output_matrix:
            
            if self.configuration['cgDNApath'] is not None:
                octave_script.write(
                    "addpath {}\n".format(self.configuration['cgDNApath']))
                
            octave_script.write("""\
sequence = '{sequence}';
params = load('{paramset}');
[nondimshapes, stiff] = constructSeqParms(sequence, params);  
save -ascii {output} stiff;""".format(
                sequence = sequence,
                paramset = self.configuration['cgDNAparamset'],
                output = output_matrix.name))
            octave_script.flush()
            
            predict  = Popen(['octave',
                                  '--no-gui',
                                  octave_script.name],
                                  stdout=PIPE, stderr=PIPE,
                                  shell=False)
            predict.wait()
            err = predict.stderr.read()
            if self._is_error(err):
                with open(octave_script.name) as error_script:
                    sys.stderr.write("SCRIPT:\n"+error_script.read())
                err = "Octave Error! Standard error stream content starts on next line:\n"+err
                raise Exception(err)
            return np.loadtxt(output_matrix.name)

    @task(sequence = IN, stiffness_matrix = OUT)
    @constraint(AppSoftware="octave,numpy,cgDNA")
    def run(self, sequence):
        """
        Return the full-oligomer stiffness matrix built using cgDNA in octave.
        """
        stiffness_matrix = self._get_stiffness_matrix(sequence)
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

    print DNAFlexRB().run(args[0])

if __name__ == "__main__":
    main()    
    
