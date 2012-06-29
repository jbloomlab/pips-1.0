"""Module for running CUPSAT ddG prediction webserver.

Written by Jesse Bloom, 2007."""
#-------------------------------------------------------------------
import os, time, urllib2, ClientForm, re, sys
#-------------------------------------------------------------------
def DDGsAlignedToSeq(protseq, pdb_matches):
    """Computes CUPSAT ddG values for a protein aligned to a PDB structure.

    This function is designed to compute ddG values for mutations to a protein
        sequence based on the alignment of that protein to a PDB file.  The
        basic idea is to use CUPSAT to compute ddG values for mutations to
        the proteins in the PDB file.  These ddG values are then assigned
        to a specified protein sequence.  If the protein sequence exactly
        matches the protein in the PDB file, then you can easily do this
        just using the 'RunCUPSAT' function.  However, this method allows
        ddG assignments to be made even if the protein sequence is not
        an exact match to the proteins in the PDB file.  In that case,
        you simply align the protein sequence with the PDB file sequence,
        and then assign the CUPSAT ddG values from the PDB file sequence
        to the corresponding aligned residues.  The actual ddG values
        are obtained by 'RunCUPSAT' using the default values for 'webserver'
        and 'print_progress'.
    'protseq' is the protein sequence for which we want to compute the ddG
        values.  It is specified as a string of one-letter amino acid codes.
    'pdb_matches' contains information about one or more PDB file chains that
        can be aligned with some or all of 'protseq'.  'pdb_matches' is a list
        of 4-tuples.  Each tuple has the following format:
        '(pdbfile, pdbchain, alignment, mapping)'.  In this tuple:
            'pdbfile' specifies the name of a PDB file (example, "1RVZ.pdb").
            'pdbchain' specifies the name of a chain in 'pdbfile' (example, "A").
                If there is only one chain in the PDB file, then 'pdbchain' should
                be 'None'.
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
    The results are returned in the dictionary 'ddgs'.  This dictionary has
        the same format as the returned variable from 'RunCUPSAT'.  Essentially,
        the dictionary is keyed by all residue numbers corresponding to positions
        in 'protseq' that correspond to a mapping in 'pdb_matches'. 'ddgs[i]'
        is the 2-tuple (wt, iddgs) where wt is the wildtype amino acid (wt = protseq[i])
        and iddgs[mut] gives the ddG for mutating from wt to some other amino acid
        mut.  If the residue wt matches that in the PDB file, then this ddG value
        is the one returned by CUPSAT for that mutation.  If the residue wt is
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
        pdb_ddgs = RunCUPSAT(pdbfile, pdbseq, pdbchain, mapping = mapping)
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
#-------------------------------------------------------------------
def RunCUPSAT(pdbfile, protseq, chain, mapping=None, webserver='http://cupsat.tu-bs.de/cupsat/custompdb.htm', print_progress=True):
    """Runs CUPSAT webserver to get ddG values for all single mutations to a protein.

    Note that this script effectively ran the CUPSAT webserver at the time that the
        script was written.  It might be rendered obsolete by modifications to
        the webserver.
    Note that these are ddG values for free energy of folding, so positive values
        are destabilizing.  CUPSAT actually returns ddG values for free energies
        of unfolding, so all of the results are multiplied by negative one.
    'pdbfile' is the path to a current PDB file for the protein in question.
    'protseq' is a string giving the sequence of the protein as one letter
        amino acid codes.  There should be no gaps or unknown residues.
    'chain' is a string giving the chain code for the protein chain in pdbfile
        for which we would like to calculate the ddG values.  If there is only
        one chain in the PDB file, then it should be set to 'None'.
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
    'webserver' specifies the web address of the CUPSAT webserver.
    'print_progress' is a Boolean switch that specifies whether we print the
        progress of the calculation.  Is 'False' by default.
    The results are returned in the dictionary 'ddgs'.  This dictionary is keyed
        by residue numbers, corresponding to the positions in 'protseq' (i.e.
        it is keyed by i where 0 <= i < len(protseq)).  'ddgs[i]' is the 
        2-tuple (wt, iddgs) where wt is the wildtype amino acid (wt = protseq[i])
        and iddgs[mut] gives the ddg for mutating from wt to some other amino acid
        mut.  If 'mapping' is specified, then the dictionary has entries only for
        residues that appear as keys in 'mapping'.
    """
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
    ddgs = {}
    for (i, pdbi) in mapping.iteritems():
        wt = protseq[i] # wildtype residue at this position
        if isinstance(pdbi, int):
            pdbi = str(pdbi)
        if print_progress:
            print "Attempting to compute ddG for PDB residue %s (%s), which corresponds to sequence residue %d of %d, of %s, chain %s..." % (pdbi, wt, i + 1, len(protseq), pdbfile, chain)
            sys.stdout.flush()
        ddgs[i] = (wt, {})
        # upload the PDB file to the webserver
        response = urllib2.urlopen(webserver)
        form = ClientForm.ParseResponse(response)[0]
        response.close()
        form.add_file(open(pdbfile), 'text/plain', pdbfile)
        request = form.click()
        # get the response after uploading the PDB file, and confirm to proceed
        response2 = urllib2.urlopen(request)
        form2 = ClientForm.ParseResponse(response2)[0]
        response2.close()
        request2 = form2.click()
        # get the response, and enter the amino acid mutation
        response3 = urllib2.urlopen(request2)
        try:
            form3 = ClientForm.ParseResponse(response3)
            form3 = form3[0]
        except:
            sys.stderr.write("form3 =\n%s" % str(form3))
            raise
        response3.close()
        form3['resno'] = pdbi # assign residue number
        form3['aawt'] = [aa_mapping[wt]] # assign residue identity
        request3 = form3.click()
        if chain:
            # get the response, and enter the chain 
            response4 = urllib2.urlopen(request3)
            form4 = ClientForm.ParseResponse(response4)[0]
            response4.close()
            form4['chainid'] = [chain]
            request4 = form4.click()
        else:
            # no need to enter chains
            response4 = urllib2.urlopen(request3)
            form4 = ClientForm.ParseResponse(response4)[0]
            response4.close()
            request4 = form4.click()
        # get the response, proceed to calculate the ddG
        response5 = urllib2.urlopen(request4)
        form5 = ClientForm.ParseResponse(response5)[0]
        response5.close()
        request5 = form5.click()
        # read the results
        response6 = urllib2.urlopen(request5)
        results = response6.read()
        response6.close()
        for mut in aa_mapping.itervalues():
            if mut == aa_mapping[wt]:
                continue # there should be no entry for the wildtype amino acid
            mutmatch = re.compile(mut + '\</td\>\n\s+\<td width\="25\%" bgcolor\="\#D9E0D8" class\="(stab|destab|middle|normal)"\> \&nbsp\;(No change|Stabilising|Destabilising)\</td\>(\n\s*)+\<td width\="25\%" bgcolor\="\#D9E0D8" class\="(stab|destab|middle|normal)"\> \&nbsp\;(Favourable|Unfavourable|No change)\</td\>(\n\s*)+\<td width\="25\%" bgcolor\="\#D9E0D8" class\="middle"\> \&nbsp\;(?P<ddg>\-{0,1}\d+\.\d{1,2})') # regular expression object to match returned ddG value
            m = mutmatch.search(results)
            if not m:
                raise ValueError, "Failed to find ddG for mutating %d (PDB %s) to %s, results are:\n\n%s" % (i, pdbi, mut, results)
            ddg = -float(m.group('ddg')) # make negative to convert to ddG for free energy of folding
            ddgs[i][1][aa_reverse_mapping[mut]] = ddg
        assert len(ddgs[i][1]) == 19
        time.sleep(0.5) # pause a bit so that we don't crash the web server
    return ddgs
#-------------------------------------------------------------------
