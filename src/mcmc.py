"""Module for performing Markov Chain Monte Carlo (MCMC) methods.

This module is written to contain the general implementation of MCMC
    methods.  The module is written as generally as possible, to allow
    it be applied to a range of different problems.

Written by Jesse Bloom, 2007."""
#-------------------------------------------------------------------
import shelve, sys, cPickle, random, os, types, copy, math
#-------------------------------------------------------------------
def MCMC(parameter_set, LogLikelihood, LogPrior, Move, nsteps, LogData, CloseLogData, last_state, output = sys.stdout, output_freq = 10, output_window = 1000, output_flush_freq = 1000, heating = False, restart = False):
    """A method for running Markov Chain Monte Carlo.

    Currently this method can perform Markov Chain Monte Carlo using the Metropolis
        acceptance criterion, but not the Metropolis-Hastings acceptance criterion.
        That is, the move proposals must be symmetric, with the chance of moving from
        x to y being the same as moving from y to x.  The method does provide options
        for performing Metropolis-coupled or "heated" MCMC.  
    The basic idea is that we are given some data, and we would
        like to construct the posterior distribution over the parameters 'parameter_set'
        in some model.  We do this using the method 'LogLikelihood' which calculates 
        Pr(data | parameter_set), the method 'LogPrior' which calculates Pr(parameter_set),
        and the method 'Move' which generates a copy of 'parameter_set' in which
        the parameter(s) have been perturbed.  Note that this perturbations must currently
        be symmetric since we have not implemented Metropolis-Hastings acceptance criterion
        for nonsymmetric moves.
    There is no returned variable from this function; instead the output is stored using
        the 'LogData' and 'CloseLogData' functions.
    'parameter_set' is a variable specifying the set of parameters for which we
        are trying to compute the posterior.  The only restrictions on this variable
        is that it must be passable to 'LogLikelihood', 'LogPrior', and 'Move'.
        Normally, 'parameter_set' is taking to specify the initial values of the
        parameters that we use to start the initial chain.  This set of parameters
        is used to start both the cold chain and all heated chains.  The only case
        when this is not the case is when 'restart' is True, in which case the value
        of 'parameter_set' is irrelevant. 
    'LogLikelihood' is a method.  'LogLikelihood(parameter_set) should return the natural log
        of the likelihood of the data given 'parameter_set'.  The data is implicitly assumed
        to be known by the function evaluating the log likelihood.
    'LogPrior' is a method.  'LogPrior(parameter_set)' should return the natural logarithm
        of the prior probability of 'parameter_set'.
    'Move' is a method that makes moves.  'Move(parameter_set)' returns a new parameter
        set 'moved_parameter_set' in which the parameters have been moved.  Note that
        'parameter_set' is unaltered in the construction of 'moved_parameter_set'.  It is
        up to you to make sure that 'Move' makes parameter moves of the appropriate size.  As
        a rule of thumb, you might want about 50-60% of the moves to be accepted.  Also,
        these moves must be symmetric, so that the probability that 'Move(parameter_set)'
        generates 'moved_parameter_set' is equal to the probability that 'Move(moved_parameter_set)'
        generates 'parameter_set'.  This symmetry is necessary to make the Metropolis acceptance
        criterion satisfy detailed balance.  
    'nsteps' is an integer specifying for how many steps the MCMC is run.  Steps run from
        0 (this is the original value of 'parameter_set' to 'nsteps').
    'LogData' is a function that is called after every step to log the current value of the
        parameters.  If there are heated chains, only the cold chain is logged.  
        'LogData' must take exactly two arguments: the first is the number of the current step
        that is being logged, and the second is the current value of the parameter set in the
        cold chain (in general, this will be an object of the same type as 'parameter_set').
        You should construct the function to somehow effectively store results from the MCMC
        chain.
    'CloseLogData' is called at the conclusion of the MCMC run.  It will be called even 
        if the run ends after an exception is called.  In general, you will want to set it
        to be a method finalize or store the data that has been logged using 'LogData'.
        It should take not arguments -- it is just called as 'CloseLogData()'.
    'last_state' is a string giving the name of a cPickle file that specifies the state of
        the MCMC chain after the last completed step.  Unpickling this object should yield
        a tuple with the first element being an integer giving the number of this step, the second
        element corresponding to the value of 'heating', the third element being the current value
        of the parameter set for the cold chain, and each remaining elements corresponding the 
        parameter set of each heated chain going from coldest to hottest.  'last_state' should
        be written even if the program crashes for some reason.  That allows the chain to be
        restarted using the 'restart' option.  The pickling is performed using the pickle
        protocol 2.
    'output' is a writable file-like object to which output is written.  By default, its value
        is standard output ('sys.stdout'), but you can set it to some other value.  If you don't
        want to write any output, set it to 'None'.  At a frequency specified by 'output_freq',
        the current posterior probability, likelihood, prior probability, acceptance frequency, and
        heated chain swapping acceptance frequency (if applicable) are written to 'output'.  The
        windows for computing the acceptance and chain-swapping acceptance are computed using a
        window of size 'output_window'.  Some preliminary and final information about the MCMC
        run are also written to 'output'.
    'output_freq' specifies the frequency with which output is written and to 'output'.
        By default, it is 10, meaning that information is written for every 10th step.  In general,
        information is only written for steps that are a multiple of 'output_freq'.
    'output_window' is used for computing the acceptance and chain swapping acceptance that 
        are written to 'output'.  These numbers are the average over the last 'output_window'
        steps.  The only exception is if the chain has just been started or restarted, in which
        case the value is 'NA' until 'output_window' steps have been performed.  The default
        value is 1000.
    'output_flush_freq' specifies the frequency with which everything written to 'output'
        is flushed.  By default, this happens every 1000 steps.
    'heating' is a variable that specifies details about the heated chains.  This is used for
        Metropolis coupling MCMC.  If you do not want to use any heated chains, then set this
        variable to 'False'.  Otherwise, it should be the 2-tuple '(nheated, heating_temp)'.
        'nheated' must be >= 1, and specifies the number of heated chains (chains in addition to
        the always present cold chain).  'heating_temp' specifies the degree of heating of these
        chains.  A heated chain is one in which the acceptance ratio (ratio of posterior probabilities)
        is raised to some power B where B < 1.  This increases the acceptance ratio for the chain.
        The heating for heated chain i (i = 1, 2,...)  is given by B_i = 1 / (1 + heating_temp * i). 
        For example, a value of 'heating' might be '(2, 0.1)'.  That means in addition to the cold chain, 
        there are heated chains with acceptance ratios raised to the powers 0.91 and 0.83.
        For heated chains, a single swap between one pair of chains is attempted at every step.
        The default value of 'heating' is 'False', meaning that no heating is done.
    'restart' is an option that can be used to restart the MCMC.  By default, this variable
        is set to 'False' which means that the MCMC is not restarted but instead is started
        fresh with the value specified by 'parameter_set'.  If you want to restart the chain,
        then set 'restart' to be the string specifying the name of the cPickle file that
        was assigned to the variable 'last_state' of the chain you are restarting.  The 
        MCMC will restart from the state specified by that file.  The value of 'heating' stored in
        'restart' must match the value of 'heated' that is used in the call.  In this case, the MCMC is 
        NOT continued for 'nsteps' new steps, but is instead continued until the restarted chain has reached 
        'nsteps' steps.  If the restarted chain has already reached or exceed 'nsteps' steps,
        then the function simply returns without doing anything. 
    """
    pickle_protocol = 2 # protocol used for all pickling
    assert isinstance(nsteps, int) and nsteps >= 1
    assert isinstance(output, file) or output == None
    assert isinstance(output_freq, int) and output_freq > 0
    assert isinstance(output_window, int) and output_window > 0
    if heating != False:
        try:
            (nheated, heating_temp) = heating
            if not (isinstance(nheated, int) and nheated >= 1):
                raise ValueError, "Invalid value of 'nheated': %s" % (str(nheated))
            if not (isinstance(heating_temp, (int, float)) and heating_temp > 0):
                raise ValueError, "Invalid value of 'heating_temp': %s" % (str(heating_temp))
        except TypeError, ValueError:
            raise ValueError, "Invalid tuple specified for 'heating': %s" % (str(heating))
    # set up the initial parameter values for the MCMC
    if restart:
        if not os.path.isfile(restart):
            raise ValueError, "Cannot find the 'restart' file of %s." % restart
        f = open(restart)
        restart_tuple = cPickle.load(f)
        f.close()
        if heating:
            if len(restart_tuple) != 3 + nheated:
                raise ValueError, "Invalid 'restart' since it does not have %d heated chains." % nheated
            (restart_step, restart_heating, restart_parameters) = (restart_tuple[0], restart_tuple[1], restart_tuple[2])
            heated_chains = [restart_tuple[3 + i] for i in range(nheated)]
        else:
            try:
                (restart_step, restart_heating, restart_parameters) = restart_tuple
            except TypeError, ValueError:
                raise ValueError, "Invalid 'restart' for no heating of %s" % (str(restart_tuple))
        if heating != restart_heating:
            raise ValueError, "'heating' is not equal to the heating specified in 'restart'."
        if restart_step >= nsteps:
            return # simply return, as we have already completed all of the steps
    if output:
        output.write("# Beginning an MCMC run.\n")
        if restart:
            output.write("# This run is a continuation of a previous run, and is being restarted from %s.\n" % restart)
            output.write("# The previous run ran for %d steps; this run will continue up until %d steps.\n" % (restart_step, nsteps))
        else:
            output.write("# This is a new run that will proceed for %d steps.\n" % (nsteps))
        output.write("# The last state of the program will be written to %s.\n" % last_state)
        if heating:
            output.write("# There will be %d heated chains, with a heating temperature of %f.\n" % (nheated, heating_temp))
        else:
            output.write("# There will be no heated chains.\n")
        output.write("# In the log that follows, updates will be printed every %d steps.\n" % output_freq)
        output.write("# Acceptance ratios will be averaged over the last %d steps.\n" % output_window)
        output.write("#\n# Definition of labels for the log are below.  Entries of NA mean that the indicated value could not be computed.\n")
        output.write("# STEP - step number.\n")
        output.write("# LOG_POSTERIOR - natural logarithm of the posterior probability of the cold chain.\n")
        output.write("# LOG_LIKELIHOOD - natural logarithm of the likelihood of the cold chain.\n")
        output.write("# LOG_PRIOR - natural logarithm of the prior of the cold chain.\n")
        output.write("# MOVED - value of 1 if the cold chain move was accepted at this step, 0 otherwise.\n")
        output.write("# SWAPPED - value of 1 if a pair of chains was swapped at this step, 0 otherwise.\n")
        output.write("# ACCEPTANCE - the acceptance frequency of moves for the cold chain averaged over the previous %d steps.\n" % output_window)
        output.write("# SWAP_ACCEPTANCE - the acceptance frequency of heated chain swaps averaged over the previous %d steps.\n" % output_window)
        output.write("#\n# STEP\tLOG_POSTERIOR\tLOG_PRIOR\tLOG_LIKELIHOOD\tMOVED\tSWAPPED\tACCEPTANCE\tSWAP_ACCEPTANCE\n")
    try: # large try... finally... loop to ensure that we close the database and write last_state
        try:
            if restart:
                istep = restart_step
                iparameters = restart_parameters
            else:
                istep = 0
                iparameters = parameter_set
                LogData(istep, iparameters) # log the parameters
            (ilogprior, iloglikelihood) = (LogPrior(iparameters), LogLikelihood(iparameters))
            ilogposterior = ilogprior + iloglikelihood
            iheated_chains = []
            if heating:
                betas = [1.0 / (1 + heating_temp * (i + 1)) for i in range(nheated)] # inverse temperatures
                if restart:
                    iheated_chains = heated_chains
                else:
                    iheated_chains = [copy.deepcopy(parameter_set) for i in range(nheated)]
                iheatedlogpriors = [LogPrior(iheatedchain) for iheatedchain in iheated_chains]
                iheatedloglikelihood = [LogLikelihood(iheatedchain) for iheatedchain in iheated_chains]
                iheatedlogposteriors = [p + l for (p, l) in zip(iheatedlogpriors, iheatedloglikelihood)]
            last_state_tuple = tuple([istep, heating, iparameters] + iheated_chains)
            if output:
                output.write("%d\t%f\t%f\t%f\tNA\tNA\tNA\tNA\n" % (istep, ilogposterior, ilogprior, iloglikelihood))
                output.flush()
            moved_list = []
            swapped_list = []
            while istep < nsteps: # main loop of the program
                # generate the moves for each chain
                nextparameters = Move(iparameters)          
                (nextlogprior, nextloglikelihood) = (LogPrior(nextparameters), LogLikelihood(nextparameters))
                nextlogposterior = nextlogprior + nextloglikelihood
                if heating:
                    nextheated_chains = []
                    nextheatedlogpriors = []
                    nextheatedloglikelihood = []
                    for iheated in range(nheated):
                        iheatedchain = iheated_chains[iheated]
                        nextheated = Move(iheatedchain)
                        nextheated_chains.append(nextheated)
                        nextheatedlogpriors.append(LogPrior(nextheated))
                        nextheatedloglikelihood.append(LogLikelihood(nextheated))
                    nextheatedlogposteriors = [p + l for (p, l) in zip(nextheatedlogpriors, nextheatedloglikelihood)]
                # decide whether to accept these moves
                r = min(1.0, math.exp(nextlogposterior - ilogposterior))
                assert 0 <= r <= 1
                if random.random() < r:
                    moved = 1
                    (iparameters, ilogprior, iloglikelihood, ilogposterior) = (nextparameters, nextlogprior, nextloglikelihood, nextlogposterior)
                else:
                    moved = 0
                if heating:
                    for iheated in range(nheated):
                        r = min(1.0, math.exp(betas[iheated] * (nextheatedlogposteriors[iheated] - iheatedlogposteriors[iheated])))
                        assert 0 <= r <= 1
                        if random.random() < r:
                            (iheated_chains[iheated], iheatedlogpriors[iheated], iheatedloglikelihood[iheated], iheatedlogposteriors[iheated]) = (nextheated_chains[iheated], nextheatedlogpriors[iheated], nextheatedloglikelihood[iheated], nextheatedlogposteriors[iheated])
                    # decide whether to swap a pair of chains
                    list_of_chains = ['cold'] + [iheated for iheated in range(nheated)]
                    random.shuffle(list_of_chains)
                    [swapj, swapk] = list_of_chains[ : 2]
                    if swapj == 'cold':
                        (logposteriorj, betaj) = (ilogposterior, 1.0)
                    else:
                        (logposteriorj, betaj) = (iheatedlogposteriors[swapj], betas[swapj])
                    if swapk == 'cold':
                        (logposteriork, betak) = (ilogposterior, 1.0)
                    else:
                        (logposteriork, betak) = (iheatedlogposteriors[swapk], betas[swapk])
                    r = min(1.0, math.exp(betaj * logposteriork + betak * logposteriorj - betaj * logposteriorj - betak * logposteriork))
                    if random.random() < r:
                        swapped = 1
                        if swapj == 'cold':
                            (tempparameters, templogposterior, temploglikelihood, templogprior) = (iparameters, ilogposterior, iloglikelihood, ilogprior)
                            (iparameters, ilogposterior, iloglikelihood, ilogprior) = (iheated_chains[swapk], iheatedlogposteriors[swapk], iheatedloglikelihood[swapk], iheatedlogpriors[swapk])
                            (iheated_chains[swapk], iheatedlogposteriors[swapk], iheatedloglikelihood[swapk], iheatedlogpriors[swapk]) = (tempparameters, templogposterior, temploglikelihood, templogprior)
                        elif swapk == 'cold':
                            (tempparameters, templogposterior, temploglikelihood, templogprior) = (iparameters, ilogposterior, iloglikelihood, ilogprior)
                            (iparameters, ilogposterior, iloglikelihood, ilogprior) = (iheated_chains[swapj], iheatedlogposteriors[swapj], iheatedloglikelihood[swapj], iheatedlogpriors[swapj])
                            (iheated_chains[swapj], iheatedlogposteriors[swapj], iheatedloglikelihood[swapj], iheatedlogpriors[swapj]) = (tempparameters, templogposterior, temploglikelihood, templogprior)
                        else:
                            (tempparameters, templogposterior, temploglikelihood, templogprior) = (iheated_chains[swapj], iheatedlogposteriors[swapj], iheatedloglikelihood[swapj], iheatedlogpriors[swapj])
                            (iheated_chains[swapj], iheatedlogposteriors[swapj], iheatedloglikelihood[swapj], iheatedlogpriors[swapj]) = (iheated_chains[swapk], iheatedlogposteriors[swapk], iheatedloglikelihood[swapk], iheatedlogpriors[swapk])
                            (iheated_chains[swapk], iheatedlogposteriors[swapk], iheatedloglikelihood[swapk], iheatedlogpriors[swapk]) = (tempparameters, templogposterior, temploglikelihood, templogprior)
                    else:
                        swapped = 0
                    swapped_list.append(swapped)
                moved_list.append(moved)
                istep += 1
                try:
                    LogData(istep, iparameters) # log the parameters
                except:
                    istep -= 1 # do this if we fail while logging the data
                    raise
                last_state_tuple = tuple([istep, heating, iparameters] + iheated_chains)
                if output and (istep % output_freq == 0):
                    if len(moved_list) >= output_window:
                        moved_list = moved_list[-output_window : ]
                        acceptance = moved_list.count(1) / float(output_window)
                        if heating:
                            swapped_list = swapped_list[-output_window : ]
                            swap_freq = swapped_list.count(1) / float(output_window)
                            output.write("%d\t%f\t%f\t%f\t%d\t%d\t%f\t%f\n" % (istep, ilogposterior, ilogprior, iloglikelihood, moved, swapped, acceptance, swap_freq))
                        else:
                            output.write("%d\t%f\t%f\t%f\t%d\tNA\t%f\tNA\n" % (istep, ilogposterior, ilogprior, iloglikelihood, moved, acceptance))
                        if not (istep % output_flush_freq):
                            output.flush()
                    else:
                        if heating:
                            output.write("%d\t%f\t%f\t%f\t%d\t%d\tNA\tNA\n" % (istep, ilogposterior, ilogprior, iloglikelihood, moved, swapped))
                        else:
                            output.write("%d\t%f\t%f\t%f\t%d\tNA\tNA\tNA\n" % (istep, ilogposterior, ilogprior, iloglikelihood, moved))
                        if not (istep % output_flush_freq):
                            output.flush()
        except:
            sys.stderr.write("MCMC ending in error after %d steps.  The final state was written to %s.\n" % (istep, last_state))
            output.write("# MCMC ending in error after %d steps.\n" % istep)
            output.write("# The final state was written to %s.\n" % (last_state))
            raise
    finally: # write last_state and close the logged data
        CloseLogData() # finalize the logged data
        try:            
            cPickle.dump(last_state_tuple, open(last_state, 'w'), protocol = pickle_protocol)
        except NameError:
            sys.stderr.write("Error, there is no final state to write to 'last_state'.\n")
    output.write("# MCMC ending normally after completing %d steps.\n" % nsteps)
    output.write("# The final state was written to %s.\n" % (last_state))
#-------------------------------------------------------------------
