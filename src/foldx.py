"""Module for running the command line version of FoldX.

This module is known to work with FoldX version 3.0 beta3 for Mac OS X,
downloaded from http://foldx.crg.es on 3-24-09.

Written by Jesse Bloom, 2009."""


import os
import re
import sys
import shutil
import tempfile
import time
from pips import ddg_inference


def DDGsAlignedToSeq(protseq, pdb_matches, foldxpath, foldxprog='FoldX.mac', hashed_ddgs=None):
    """Computes FoldX ddG values for a protein aligned to a PDB structure.

    This function is designed to compute ddG values for mutations to a protein
        sequence based on the alignment of that protein to a PDB file.  The
        basic idea is to use FoldX to compute ddG values for mutations to
        the proteins in the PDB file.  These ddG values are then assigned
        to a specified protein sequence.  If the protein sequence exactly
        matches the protein in the PDB file, then you can easily do this
        just using the 'RunFoldX' function.  However, this method allows
        ddG assignments to be made even if the protein sequence is not
        an exact match to the proteins in the PDB file.  In that case,
        you simply align the protein sequence with the PDB file sequence,
        and then assign the FoldX ddG values from the PDB file sequence
        to the corresponding aligned residues.  The actual ddG values
        are obtained by 'RunFoldX'.
    'protseq' is the protein sequence for which we want to compute the ddG
        values.  It is specified as a string of one-letter amino acid codes.
    'pdb_matches' contains information about one or more PDB file chains that
        can be aligned with some or all of 'protseq'.  'pdb_matches' is a list
        of 4-tuples.  Each tuple has the following format:
        '(pdbfile, pdbchain, alignment, mapping)'.  In this tuple:
            'pdbfile' specifies the name of a PDB file (example, "1RVZ.pdb").
            'pdbchain' specifies the name of a chain in 'pdbfile' (example, "A").
            'alignment' specifies an alignment of 'protseq' with the sequence of
                the chain 'pdbchain' in 'pdbfile'.  This alignment should be in the
                form of a 2-tuple '(alignedprotseq, alignedpdbseq)'.  This 2-tuple
                represents 'protseq' and the sequence of 'pdbchain' aligned to
                each other, with gap ('-') characters inserted as appropriate
                to make the sequence align.  Stripping all of the gaps out of
                'alignedprotseq' should yield 'protseq'.  Stripping all of the
                gaps out of 'alignedpdbseq' should yield the sequence of
                the chain 'pdbchain' in 'pdbfile'.  Call this latter stripped
                sequence 'pdbseq'.  Crucially, 'pdbseq' must exactly correspond
                to the sequence in the PDB file.
            'mapping' maps the residues in 'pdbseq' to the corresponding residue
                numbers in the PDB file.  Specifically, for each residue i
                in 'pdbseq' (0 <= i < len(pdbseq)), there can be a key in 
                'mapping'.  The value corresponding to this key should be a
                string giving the number of the residue in that chain 'pdbchain'
                that corresponds to this residue.  Any residues in 'pdbseq'
                that do not have keys in 'mapping' do not have ddG values computed.
                If 'mapping' is equal to 'None', then residue i is mapped to the
                number i + 1 for all i.
        'pdb_matches' can contain multiple 4-tuples if different parts of 'protseq' 
            are to be matched to different PDB files or different chains in the
            same PDB file.  However, in this case each residue in 'protseq' can
            only be aligned to one residue; it cannot be aligned to residues in
            more than one 4-tuple entry of 'pdb_matches'.
    'foldxpath' specifies the path to the FoldX executable, which has a name
        given by 'foldxprog'.  This directory must also contain the file "rotabase.txt".
    'foldxprog' gives the name of the FoldX executable.
    'hashed_ddgs' is an optional switch that allows us to store ddG values from
        partiall complete calculations.  By default, it is 'None', meaning that 
        no hashing is done.  It can be set to a string specifying a valid 
        filename.  In this case, if the program terminates before all ddG
        values have been computed, the values computed so far are written
        to the file hashed_ddgs.  If the function is called again with
        the same value of hashed_ddgs and the file of that name can be 
        found, then those already computed values are used, and just 
        the remaining non-hashed ddG values are computed.  When the program
        terminates normally, any existing hashed_ddgs file is deleted.
    The results are returned in the dictionary 'ddgs'.  This dictionary has
        the same format as the returned variable from 'RunFoldX'.  Essentially,
        the dictionary is keyed by all residue numbers corresponding to positions
        in 'protseq' that correspond to a mapping in 'pdb_matches'. 'ddgs[i]'
        is the 2-tuple (wt, iddgs) where wt is the wildtype amino acid (wt = protseq[i])
        and iddgs[mut] gives the ddG for mutating from wt to some other amino acid
        mut.  If the residue wt matches that in the PDB file, then this ddG value
        is the one returned by FoldX for that mutation.  If the residue wt is
        different from that in the PDB file, then the ddG value is calculated
        by adding the ddG value for mutating the PDB residue (wt_pdb) to mut
        to the negative of the value of mutating wt_pdb to wt.  
    """
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    assert isinstance(protseq, str)
    ddgs = {}
    for (pdbfile, pdbchain, (alignedprotseq, alignedpdbseq), mapping) in pdb_matches:
        assert os.path.isfile(pdbfile)
        assert len(alignedprotseq) == len(alignedpdbseq)
        assert protseq == alignedprotseq.replace('-', '')
        pdbseq = alignedpdbseq.replace('-', '')
        pdb_ddgs = RunFoldX(pdbfile, pdbseq, pdbchain, foldxpath, mapping=mapping, foldxprog=foldxprog, hashed_ddgs=hashed_ddgs)
        iprot = ipdb = 0
        for (protres, pdbres) in zip(alignedprotseq, alignedpdbseq):
            if (protres != '-') and (pdbres != '-') and ipdb in mapping:
                if iprot in ddgs:
                    raise ValueError, "Duplicate ddG values for residue %d in protseq." % iprot
                if ipdb not in pdb_ddgs:
                    raise ValueError, "No key in PDB ddGs for ipdb = %d (iprot = %d)." % (ipdb, iprot)
                if protres == pdbres: # residues are identical in protseq and PDB file
                    ddgs[iprot] = pdb_ddgs[ipdb]
                else: # residues differ in protseq and PDB file
                    iddgs = {}
                    protres_to_pdbres = -pdb_ddgs[ipdb][1][protres]
                    for mut in amino_acids:
                        if mut == protres:
                            continue
                        elif mut == pdbres:
                            pdbres_to_mut = 0.0
                        else:
                            pdbres_to_mut = pdb_ddgs[ipdb][1][mut]
                        iddgs[mut] = protres_to_pdbres + pdbres_to_mut
                    ddgs[iprot] = (protres, iddgs)
                assert protseq[iprot] == protres == ddgs[iprot][0]
            if protres != '-':
                iprot += 1
            if pdbres != '-':
                ipdb += 1
    return ddgs



def RunFoldX(pdbfile, protseq, chain, foldxpath, mapping=None, foldxprog='FoldX.mac', print_progress=True, repair_pdb=True, ionic_strength=0.05, hashed_ddgs=None):
    """Runs FoldX to get ddG values for all single mutations to a protein.

    This function uses the command line version of FoldX to compute the ddG values
        for single mutations to a protein.  It does this with the position scan feature.
        It is tested with FoldX version 3.0 beta for Mac OS X.
    Note that these are ddG values for free energy of folding, so positive values
        are destabilizing.  
    'pdbfile' is the path to a current PDB file for the protein in question.
    'protseq' is a string giving the sequence of the protein as one letter
        amino acid codes.  There should be no gaps or unknown residues.
    'chain' is a string giving the chain code for the protein chain in pdbfile
        for which we would like to calculate the ddG values.  If 'chain' has a 
        value of 'None', it is automatically changed to "A".
    'foldxpath' specifies the path to the FoldX executable, which has a name
        given by 'foldxprog'.  This directory must also contain the file "rotabase.txt".
    'mapping' specifies the relationship between 'protseq' and the relevant chain in 
        'pdbfile'.  If it has its default value of 'None', then 'protseq' and 
        the chain in 'pdbfile' specify the exact same protein sequence, with
        the residue 'protseq[i]' corresponding to the residue numbered i + 1
        in 'pdbfile'.  If this mapping is not exact, then it is necessary to
        remap the residues.  In this case, 'mapping' should be a dictionary keyed
        by integers.  'mapping[i]' should be a string giving the residue number in 
        'pdbfile' that corresponds to residue 'i' in 'protseq', where 'i' can
        be any value satisfying 0 <= i < len(protseq).  Only residues that
        appear in 'mapping' then have their ddG values computed.  'mapping[i]'
        can also be an integer corresponding to the residue number.
    'foldxprog' gives the name of the FoldX executable.
    'print_progress' is a Boolean switch specifying that we print progress updates
        to standard output.
    'repair_pdb' is a Boolean switch specifying whether we run the FoldX RepairPDB
        command prior to performing the mutagenesis.  This is only done if
        'repair_pdb' has its default value of True.
    'ionic_strength' specifies the ionic strength of the solution in M.  This value
        is used by FoldX in its stability calculations.
    'hashed_ddgs' is an optional switch that allows us to store ddG values from
        partiall complete calculations.  By default, it is 'None', meaning that 
        no hashing is done.  It can be set to a string specifying a valid 
        filename.  In this case, if the program terminates before all ddG
        values have been computed, the values computed so far are written
        to the file hashed_ddgs.  If the function is called again with
        the same value of hashed_ddgs and the file of that name can be 
        found, then those already computed values are used, and just 
        the remaining non-hashed ddG values are computed.  When the program
        terminates normally, any existing hashed_ddgs file is deleted.
    The results are returned in the dictionary 'ddgs'.  This dictionary is keyed
        by residue numbers, corresponding to the positions in 'protseq' (i.e.
        it is keyed by i where 0 <= i < len(protseq)).  'ddgs[i]' is the 
        2-tuple (wt, iddgs) where wt is the wildtype amino acid (wt = protseq[i])
        and iddgs[mut] gives the ddg for mutating from wt to some other amino acid
        mut.  If 'mapping' is specified, then the dictionary has entries only for
        residues that appear as keys in 'mapping'.
    """
    assert hashed_ddgs == None or isinstance(hashed_ddgs, str)
    if chain == None:
        chain = 'A'
    aa_mapping = {'A':'ALA', 'C':'CYS', 'D':'ASP', 'E':'GLU', 'F':'PHE', 'G':'GLY', 'H':'HIS', 
                  'I':'ILE', 'K':'LYS', 'L':'LEU', 'M':'MET', 'N':'ASN', 'P':'PRO', 'Q':'GLN',
                  'R':'ARG', 'S':'SER', 'T':'THR', 'V':'VAL', 'W':'TRP', 'Y':'TYR'}
    aa_reverse_mapping = {}
    for (one, three) in aa_mapping.iteritems():
        aa_reverse_mapping[three] = one
    assert os.path.isfile(pdbfile)
    if mapping == None:
        # construct the default mapping
        mapping = {}
        for i in range(len(protseq)):
            mapping[i] = str(i + 1)
    elif not isinstance(mapping, dict):
        raise ValueError, 'mapping is not a dictionary.'
    if not os.path.isdir(foldxpath):
        raise IOError("Cannot find 'foldxpath' of %s." % foldxpath)
    if not os.path.isfile("%s/%s" % (foldxpath, foldxprog)):
        raise IOError("Cannot find 'foldxprog' of %s." % foldxprog)
    if not os.path.isfile("%s/rotabase.txt" % foldxpath):
        raise IOError("Cannot find 'rotabase.txt'.")
    ddgs = {}
    if hashed_ddgs and os.path.isfile(hashed_ddgs):
        (datetime, ddgs) = ddg_inference.ReadDDGs(hashed_ddgs)
        if print_progress:
            print "Using hashed ddG values read from %s." % hashed_ddgs
            print "These ddG values are described as %s" % datetime
            print "There are hashed ddG values for %d residues." % len(ddgs)
        os.remove(hashed_ddgs)
    # create a temporary directory for running the program
    currdir = os.getcwd()
    tempdir = tempfile.mkdtemp()
    try:
        # do stuff in temporary directory
        shutil.copy(pdbfile, tempdir)
        shutil.copy("%s/rotabase.txt" % foldxpath, tempdir)
        os.chdir(tempdir)
        if repair_pdb:
            # repair the PDB file using FoldX
            if print_progress:
                print "First repairing the PDB file with FoldX..."
                sys.stdout.flush()
            open('run.txt', 'w').write("<TITLE>FOLDX_runscript;\n<JOBSTART>#;\n<PDBS>%s;\n<BATCH>#;\n<COMMANDS>FOLDX_commandfile;\n<RepairPDB>#;\n<END>#;\n<OPTIONS>FOLDX_optionfile;\n<Temperature>298;\n<R>#;\n<pH>7;\n<IonStrength>0.050;\n<water>-CRYSTAL;\n<metal>-CRYSTAL;\n<VdWDesign>2;\n<OutPDB>false;\n<pdb_hydrogens>false;\n<complex_with_DNA> true;\n<END>#;\n<JOBEND>#;\n<ENDFILE>#;" % (pdbfile))
            os.system("%s/%s -runfile run.txt >& temp.txt" % (foldxpath, foldxprog))
            pdbfile = "RepairPDB_%s" % pdbfile
            if not os.path.isfile(pdbfile):
                raise ValueError("Cannot find repaired PDB.")
        pdb_id = os.path.splitext(pdbfile)[0]
        for (i, pdbi) in mapping.iteritems():
            wt = protseq[i] # wildtype residue at this position
            if isinstance(pdbi, int):
                pdbi = str(pdbi)
            if i in ddgs and len(ddgs[i][1]) == 19:
                if print_progress:
                    print "ddGs for PDB residue %s (%s), which corresponds to sequence residue %d of %d, of %s, chain %s, have already been computed and are being read from hashed values." % (pdbi, wt, i + 1, len(protseq), pdbfile, chain)
                if wt != ddgs[i][0]:
                    raise ValueError("Residue identity mismatch with hashed ddGs for residue %d (%s versus %s)" % (i, wt, ddgs[i][0]))
                continue
            if print_progress:
                print "Attempting to compute ddGs for PDB residue %s (%s), which corresponds to sequence residue %d of %d, of %s, chain %s..." % (pdbi, wt, i + 1, len(protseq), pdbfile, chain)
                sys.stdout.flush()
            ddgs[i] = (wt, {})
            # prepare the run file
            open('run.txt', 'w').write("<TITLE>FOLDX_runscript;\n<JOBSTART>#;\n<PDBS>%s;\n<BATCH>#;\n<COMMANDS>FOLDX_commandfile;\n<PositionScan>#,%s%s%sa;\n<END>#;\n<OPTIONS>FOLDX_optionfile;\n<Temperature>298;\n<R>#;\n<pH>7;\n<IonStrength>%.3f;\n<water>-CRYSTAL;\n<metal>-CRYSTAL;\n<VdWDesign>2;\n<OutPDB>false;\n<pdb_hydrogens>false;\n<complex_with_DNA> true;\n<END>#;\n<JOBEND>#;\n<ENDFILE>#;" % (pdbfile, wt, chain, pdbi, ionic_strength))
            os.system("%s/%s -runfile run.txt >& temp.txt" % (foldxpath, foldxprog))
            resultfile = 'energies_%s_%s.txt' % (pdbi, pdb_id)
            if not os.path.isfile(resultfile):
                raise IOError("Failed to find result file of %s." % resultfile)
            results = open(resultfile).read()
            dgs = {}
            for mut in aa_mapping.itervalues():
                mutmatch = re.compile("%s_%s_%s.txt\t(?P<dg>\-{0,1}\d+\.\d+)\t" % (mut, pdbi, pdb_id))
                m = mutmatch.search(results)
                if not m:
                    raise ValueError("Failed to find dG for residue %s at position %s of %s; the results are:\n%s" % (mut, pdbi, pdb_id, results))
                dgs[aa_reverse_mapping[mut]] = float(m.group('dg')) # free energy of folding for this mutant
            for (mut, dg) in dgs.iteritems():
                if mut != wt:
                    ddg = dgs[mut] - dgs[wt]
                    ddgs[i][1][mut] = ddg
    finally:
        os.chdir(currdir) # return to current directory
        for file in os.listdir(tempdir):
            os.remove("%s/%s" % (tempdir, file))
        os.rmdir(tempdir)
        if hashed_ddgs:
            ddg_inference.WriteDDGs(ddgs, hashed_ddgs, 'FOLDX ddG values hashed at %s.' % time.asctime())
    return ddgs
