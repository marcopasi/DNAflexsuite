# File: DNAFlexibility.py 
#
# Copyright (C) 2016 Marco Pasi <marco.pasi@nottingham.ac.uk> 
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
#"DNAFlexibility.py v0.1 (C) 2016 Marco Pasi"
#

from mug import datatypes as mug_datatypes
from mugtools.resource import Resource
from mugtools.tool import Tool

from pycompss.api.task import task
from pycompss.api.constraint import constraint
from pycompss.api.parameter import IN, OUT

from itertools import imap

#------------------------------------------------------------------------------
def ogrouper(seq, size):
    """
    Yield overlapping groupings of size =size= from sequence =seq=.
    Example:
    >>>map("".join, list(ogrouper("abcdef", 3)))
    ['abc', 'bcd', 'cde', 'def']
    """
    it = iter(seq)
    stack = (None,)
    for n in xrange(size-1):
        stack += (it.next(),)        
    while True:
        values = stack[1:]
        values += (it.next(),)
        stack = values
        yield values

#------------------------------------------------------------------------------
class DNAFlexibility(Tool):
    """
    Abstract class defining a predictor of DNA flexibility as a function of
    of sequence, based on the assumption that flexibility's sequence-dependence
    can be accurately described at a given k-nucleotide level. 
    The particular form of the prediction output must be defined in subclasses. 
    """
    
    input_data_type = mug_datatypes.Sequence
    output_data_type = mug_datatypes.Undefined # must be defined in subclasses
    configuration = {}
    
    def _get_flexibility(self, kad):
        """Return the flexibility of the specified kad."""
        pass
        
    def _get_flexibilities(self, sequence, k=4):
        """Cycle through k-ads and call _get_flexibility on each."""
        for kad in imap("".join, ogrouper(sequence, k)):
            yield self._get_flexibility(kad)

    def run(self):
        pass
