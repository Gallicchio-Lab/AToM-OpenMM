# File Based Replica Exchange class
"""
The core module of ASyncRE-OpenMM: a framework to run asynchronous replica exchange calculations with OpenMM

Authors:
Emilio Gallicchio <emilio.gallicchio@gmail.com>

"""
from __future__ import print_function
from __future__ import division
import os
import sys
import time
import pickle
import random, glob
import shutil
import logging, logging.config
import signal

from atom_openmm.gibbs_sampling import *

from atom_openmm.ommreplica import *
from atom_openmm.ommworker import *
from atom_openmm.local_openmm_transport import *
from atom_openmm.utils.config import parse_config

import multiprocessing as mp

__version__ = '8.2.1'

class JobManager(object):
    """
    Class to set up and run asynchronous file-based RE calculations
    """
    logging.config.fileConfig(os.path.join(os.path.dirname(__file__), "utils/logging.conf"))

    def __init__(self, command_file, options):
        try:
            import setproctitle
            setproctitle.setproctitle("AToM %s" % command_file)
        except:
            pass

        self._setLogger()

        self.command_file = command_file
        if not os.path.exists(self.command_file):
           self._exit('No such file: %s'%self.command_file)

        self.jobname = os.path.splitext(os.path.basename(command_file))[0]
        self.keywords = parse_config(self.command_file)

        self._checkInput()
        self._printStatus()

        #set to False to run without exchanges
        self.exchange = True

        # Set the async mode to True by default
        self.async_mode = self.keywords.get('ASYNC_MODE', True)

        if self.async_mode:
            #catch ctrl-C and SIGTERM to terminate threads gracefully
            signal.signal(signal.SIGINT, self._signal_handler)
            signal.signal(signal.SIGTERM, self._signal_handler)

    def _exit(self, message):
        """Print and flush a message to stdout and then exit."""
        self.checkpointJob()
        self._cleanup()
        self.logger.info(message)
        sys.stdout.flush()
        if not self.async_mode:
            sys.exit(1)

        self.logger.info('Waiting for children to complete ...')
        while True:
            time.sleep(1)
            if not mp.active_children():
                break
        time.sleep(20)
        sys.exit(1)


    def getVersion(self):
        return __version__

    def _cleanup(self):
        try:
            for worker in self.openmm_workers:
                worker.finish()
        except:
            pass

    def _signal_handler(self, sig, frame):
        msg = "Termination signal %s detected ... cleaning up." % str(sig)
        self._exit(msg)

    def _setLogger(self):
        self.logger = logging.getLogger("async_re")

    def __getattribute__(self, name):
        if name == 'replicas_waiting':
            # Return a list of replica indices of replicas in a wait state.
            return [k for k in range(self.nreplicas)
                    if self.status[k]['running_status'] == 'W']
        elif name == 'states_waiting':
            # Return a list of state ids of replicas in a wait state.
            return [self.status[k]['stateid_current']
                    for k in self.replicas_waiting]
        elif name == 'replicas_waiting_to_exchange':
            # Return a list of replica indices of replicas in a wait state that
            # have ALSO completed at least one cycle.
            return [k for k in range(self.nreplicas)
                    if (self.status[k]['running_status'] == 'W' and
                        self.status[k]['cycle_current'] > 1)]
        elif name == 'states_waiting_to_exchange':
            # Return a list of state ids of replicas in a wait state that have
            # ALSO completed at least one cycle.
            return [self.status[k]['stateid_current']
                    for k in self.replicas_waiting_to_exchange]
        elif name == 'waiting':
            return len(self.replicas_waiting)
        elif name == 'replicas_running':
            # Return a list of replica indices of replicas in a running state.
            return [k for k in range(self.nreplicas)
                    if self.status[k]['running_status'] == 'R']
        elif name == 'running':
            return len(self.replicas_running)
        else:
            return object.__getattribute__(self,name)

    def _printStatus(self):
        """Print a report of the input parameters."""
        self.logger.info("ASyncRE-OpenMM, Version %s" % self.getVersion())
        self.logger.info("command_file = %s", self.command_file)
        self.logger.info("jobname = %s", self.jobname)
        self.logger.info("Keywords:")
        for k,v in self.keywords.items():
            self.logger.info("%s: %s", k, v)

    def _checkInput(self):
        """
        Check that required parameters are specified. Parse these and other
        optional settings.
        """
        # Required Options
        #
        # basename for the job
        self.basename = self.keywords.get('BASENAME')
        if self.basename is None:
            self._exit('BASENAME needs to be specified')

        #job transport mechanism
        self.transport_mechanism = self.keywords.get('JOB_TRANSPORT')
        if self.transport_mechanism is None:
            self._exit('JOB_TRANSPORT needs to be specified')
        #only LOCAL_OPENMM is supported for now
        if self.transport_mechanism != "LOCAL_OPENMM" :
            self._exit("unknown JOB_TRANSPORT %s" % self.transport_mechanism)
        # reset job transport
        self.transport = None

        if self.transport_mechanism == "LOCAL_OPENMM":
            if self.keywords.get('NODEFILE') is None:
                self._exit("NODEFILE needs to be specified")
            nodefile = self.keywords.get('NODEFILE')
            """
            check the information in the nodefile. there should be six columns in the  
            nodefile. They are 'node name', 'slot number', 'number of threads', 
            'platform','username', and 'name of the temperary folder'
            """
            node_info= []
            try:
                f = open(nodefile, 'r')
                line=f.readline()
                nodeid = 0
                while line:
                    lineID=line.split(",")
                    node_info.append({})
                    node_info[nodeid]["node_name"] = str(lineID[0].strip())
                    node_info[nodeid]["slot_number"] = str(lineID[1].strip())
                    node_info[nodeid]["threads_number"] = str(lineID[2].strip())
                    node_info[nodeid]["arch"] = str(lineID[3].strip())
                    node_info[nodeid]["user_name"] = str(lineID[4].strip())
                    node_info[nodeid]["tmp_folder"] = str(lineID[5].strip())
                    #tmp_folder has to be pre-assigned
                    if node_info[nodeid]["tmp_folder"] == "":
                        self._exit('tmp_folder in nodefile needs to be specified')
                    nodeid += 1
                    line = f.readline()
                f.close()
            except:
                self._exit("Unable to process nodefile %s" % nodefile)
                # reset job transport
                self.transport = None
            #set the nodes information
            self.num_nodes = len(node_info)
            self.compute_nodes=node_info
            #Can print out here to check the node information
            self.logger.info("compute nodes: %s", ', '.join([n['node_name'] for n in node_info]))

        # execution time in minutes
        self.walltime = float(self.keywords.get('WALL_TIME'))
        if self.walltime is None:
            self._exit('WALL_TIME (in minutes) needs to be specified')

        # Optional variables
        #
        env = self.keywords.get('ENGINE_ENVIRONMENT')
        if env is not None and env != '':
            self.engine_environment = env.split(',')
        else:
            self.engine_environment = []

        # number of replicas (may be determined by other means)
        self.nreplicas = None

        # verbose printing
        self.verbose = self.keywords.get('VERBOSE', False)
        if self.verbose and self.logger:
            self.logger.setLevel(logging.DEBUG)

        self.implicitsolvent =  self.keywords.get('IMPLICITSOLVENT')
        self.totalsteps = self.keywords.get('PRODUCTION_STEPS')
        self.jobname = self.keywords.get('BASENAME')
        self.stepgap = self.keywords.get('PRNT_FREQUENCY')

    def setupJob(self):
        self.transport = LocalOpenMMTransport(self.basename, self.openmm_workers, self.openmm_replicas)
        # create status table
        self.status = [{'stateid_current': k, 'running_status': 'W',
                        'cycle_current': 1} for k in range(self.nreplicas)]
        for replica in self.openmm_replicas:
            self.status[replica._id]['cycle_current'] = replica.get_cycle()
            self.status[replica._id]['stateid_current'] = replica.get_stateid()
            self.logger.info("Replica %d Cycle %d Stateid %d" % (replica._id, self.status[replica._id]['cycle_current'], self.status[replica._id]['stateid_current']))
        self.updateStatus()

        self.print_status()

    def scheduleJobs(self):
        # Gets the wall clock time for a replica to complete a cycle
        # If unspecified it is estimated as 1% of job wall clock time
        replica_run_time = self.keywords.get('REPLICA_RUN_TIME')
        if self.keywords.get('REPLICA_RUN_TIME') is None:
            replica_run_time = int(round(self.walltime/100.))
        else:
            replica_run_time = int(self.keywords.get('REPLICA_RUN_TIME'))
        # double it to give time for current running processes
        # and newly submitted processes to complete
        replica_run_time *= 2

        # Time in between cycles in seconds
        # If unspecified it is set as 30 secs
        if self.keywords.get('CYCLE_TIME') is None:
            cycle_time = 30.0
        else:
            cycle_time = float(self.keywords.get('CYCLE_TIME'))

        if self.keywords.get('MIN_TIME') is None:
            min_time = 1
        else:
            min_time = float(self.keywords.get('MIN_TIME'))

        if self.keywords.get('CHECKPOINT_TIME') is None:
            checkpoint_time = cycle_time
        else:
            checkpoint_time = float(self.keywords.get('CHECKPOINT_TIME'))

        sample_steps = int(self.keywords.get('PRNT_FREQUENCY'))
        cycle_steps = int(self.keywords.get('PRODUCTION_STEPS'))
        checkpoint_frequency = int(self.keywords.get('CHECKPOINT_FREQUENCY', 0))
        cycle_to_sample = sample_steps//cycle_steps
        enough_samples = False
        max_samples = None

        def current_samples():
            return [(replica.get_cycle()-1)//cycle_to_sample for replica in self.openmm_replicas]

        if self.keywords.get('MAX_SAMPLES')  is not None:
            max_samples_str = self.keywords.get('MAX_SAMPLES')
            max_samples = int(max_samples_str)
            starting_samples = 0

            if isinstance(max_samples_str, str) and max_samples_str.startswith("+"):
                # Handle cases where we want to increase the number of samples from a starting checkpoint
                if not os.path.isfile("starting_sample"):
                    with open("starting_sample", "w") as f:
                        f.write(f"{min(current_samples())}\n")
                with open("starting_sample", "r") as f:
                    starting_samples = int(f.read().strip())
                    max_samples += starting_samples

            self.logger.info(f"Target number of samples: {max_samples}. Current samples: {min(current_samples())}")

            enough_samples = all( [ repl_sample >= max_samples for repl_sample in current_samples() ] )
            if enough_samples:
                self.logger.info("All replicas collected the requested number of samples (%d)" % max_samples)

        start_time = time.time()
        end_time = start_time + 60*(self.walltime - replica_run_time)
        last_checkpoint_time = start_time

        while ( time.time() < end_time and
                self.transport.numNodesAlive() > 0 and
                not enough_samples ) :
            if max_samples is not None:
                with open("progress", "w") as f:
                    f.write(f"{(min(current_samples())-starting_samples) / (max_samples - starting_samples)}\n")
            current_time = time.time()

            self.updateStatus()
            self.print_status()
            self.launchJobs()
            self.updateStatus()
            self.print_status()

            self.transport.ProcessJobQueue(min_time,cycle_time)

            self.updateStatus()
            self.print_status()
            if self.exchange:
                self.doExchanges()
            self._write_status()
            self.print_status()
            self.transport.fixnodes()

            if current_time - last_checkpoint_time > checkpoint_time or (checkpoint_frequency != 0 and all([(sample % checkpoint_frequency) == 0 for sample in current_samples()] )):
                self.logger.info("Checkpointing ...")
                self.checkpointJob()
                last_checkpoint_time = current_time
                self.logger.info("done.")

            #terminates if enough samples have been collected
            if max_samples is not None:
                enough_samples = all( [ (replica.get_cycle()-1)/cycle_to_sample >=  max_samples
                                        for replica in self.openmm_replicas ]  )
                if enough_samples:
                    self.logger.info("All replicas collected the requested number of samples (%d)" % max_samples)

        if self.transport.numNodesAlive() <= 0 :
            self.logger.info("No compute devices are alive. Quitting.")
        else:
            if time.time() >= end_time :
                self.logger.info("Requested simulation time completed (%d mins)" % self.walltime)
            self.logger.info("Normal termination.")

        self.transport.DrainJobQueue()
        self.updateStatus()
        self.print_status()
        self.waitJob()
        self.checkpointJob()
        self.cleanJob()

    def waitJob(self):
        # wait until all jobs are complete
        completed = False
        while not completed:
            self.updateStatus()
            completed = True
            for k in range(self.nreplicas):
                if self.status[k]['running_status'] == "R":
                    completed = False
            time.sleep(1)

    def cleanJob(self):
        self._cleanup()
        return

    def checkpointJob(self):
        #defined in subclasses
        pass

    def _write_status(self):
        pass

    def _read_status(self):
        pass

    def print_status(self):
        """
        Writes to BASENAME_stat.txt a text version of the status of the RE job.
        It's fun to follow the progress in real time by doing:
        watch cat BASENAME_stat.txt
        """
        log = 'Replica  State  Status  Cycle \n'
        for k in range(self.nreplicas):
            log += ('%6d   %5d  %5s  %5d \n'%
                    (k,self.status[k]['stateid_current'],
                     self.status[k]['running_status'],
                     self.status[k]['cycle_current']))
        log += 'Running = %d\n'%self.running
        log += 'Waiting = %d\n'%self.waiting

        logfile = '%s_stat.txt'%self.basename
        ofile = _open(logfile,'w')
        ofile.write(log)
        ofile.close()

    def _buildInpFile(self, repl):
        pass

    def updateStatus(self):
        """Scan the replicas and update their states."""
        self.transport.poll()
        for k in range(self.nreplicas):
            self._updateStatus_replica(k)
        self._write_status()

    def _updateStatus_replica(self, replica):
        this_cycle = self.status[replica]['cycle_current']
        if self.status[replica]['running_status'] == 'R':
            if self.transport.isDone(replica,this_cycle):
                self.status[replica]['running_status'] = 'S'
                #MD engine modules implement ways to check for completion.
                #by testing existence of output file, etc.
                if self._hasCompleted(replica,this_cycle):
                    self.status[replica]['cycle_current'] += 1
                else:
                    self.logger.warning('_updateStatus_replica(): restarting replica %s (cycle %s)',
                                        replica, this_cycle)
                self.status[replica]['running_status'] = 'W'
        self.update_state_of_replica(replica)

    def _njobs_to_run(self):
        # size of subjob buffer as a percentage of job slots
        subjobs_buffer_size = self.keywords.get('SUBJOBS_BUFFER_SIZE')
        if subjobs_buffer_size is None:
            subjobs_buffer_size = 0.5
        else:
            subjobs_buffer_size = float(subjobs_buffer_size)

        # launch new replicas if the number of submitted/running subjobs is
        # less than the number of available slots
        # (total_cores/subjob_cores) + 50%
        available_slots = self.num_nodes
        max_njobs_submittable = int((1.+subjobs_buffer_size)*available_slots)
        nlaunch = self.waiting - max(2,self.nreplicas - max_njobs_submittable)
        nlaunch = max(0,nlaunch)
        if self.verbose:
            self.logger.debug('available_slots: %d', available_slots)
            self.logger.debug('max job queue size: %d', max_njobs_submittable)
            self.logger.debug('running/submitted subjobs: %d', self.running)
            self.logger.debug('waiting replicas: %d', self.waiting)
            self.logger.debug('replicas to launch: %d', nlaunch)
        return nlaunch

    def _cycle_of_replica(self,repl):
        return self.status[repl]['cycle_current']

    def launchJobs(self):
        """
        Scans the replicas in wait state and randomly launches them
        """
        if self.async_mode:
            jobs_to_launch = self._njobs_to_run()
            if jobs_to_launch > 0:
                #prioritize replicas that are most behind
                wait = sorted(self.replicas_waiting, key=self._cycle_of_replica)
                #  random.shuffle(wait)
                n = min(jobs_to_launch,len(wait))
                for k in wait[0:n]:
                    self.logger.info('Launching replica %d cycle %d', k, self.status[k]['cycle_current'])
                    # the _launchReplica function is implemented by
                    # MD engine modules
                    status = self._launchReplica(k,self.status[k]['cycle_current'])
                    if status != None:
                        self.status[k]['running_status'] = 'R'
        else:
            # In sync mode launch all replicas at once and they will be executed sequentially
            for k in range(self.nreplicas):
                self.logger.info('Launching replica %d cycle %d', k, self.status[k]['cycle_current'])
                # the _launchReplica function is implemented by
                # MD engine modules
                status = self._launchReplica(k,self.status[k]['cycle_current'])
                if status is not None:
                    self.status[k]['running_status'] = 'R'

    def doExchanges(self):
        """Perform exchanges among waiting replicas using Gibbs sampling."""

        #check the exchange
        if self.verbose:
            self.logger.info("starting the replica exchange")

        replicas_to_exchange = self.replicas_waiting_to_exchange
        states_to_exchange = self.states_waiting_to_exchange
        nreplicas_to_exchange = len(replicas_to_exchange)
        if nreplicas_to_exchange < 2:
            return 0

        if self.verbose:
            self.logger.debug('Initiating exchanges amongst %d replicas:', nreplicas_to_exchange)

        exchange_start_time = time.time()

        # Matrix of replica energies in each state.
        # The computeSwapMatrix() function is defined by application classes
        matrix_start_time = time.time()
        swap_matrix = self._computeSwapMatrix(replicas_to_exchange,
                                              states_to_exchange)
        if self.verbose:
            self.logger.info(swap_matrix)
        matrix_time = time.time() - matrix_start_time

        sampling_start_time = time.time()

        for repl_i in replicas_to_exchange:
            #repl_i = choice(replicas_to_exchange)
            sid_i = self.status[repl_i]['stateid_current']
            curr_states = [self.status[repl_j]['stateid_current']
                           for repl_j in replicas_to_exchange]
            repl_j = pairwise_independence_sampling(repl_i,sid_i,
                                                    replicas_to_exchange,
                                                    curr_states,
                                                    swap_matrix)
            if repl_j != repl_i:
                sid_i = self.status[repl_i]['stateid_current']
                sid_j = self.status[repl_j]['stateid_current']
                self.status[repl_i]['stateid_current'] = sid_j
                self.status[repl_j]['stateid_current'] = sid_i
                self.logger.info("Replica %d new state %d" % (repl_i, sid_j))
                self.logger.info("Replica %d new state %d" % (repl_j, sid_i))

        # Uncomment to debug Gibbs sampling:
        # Actual and observed populations of state permutations should match.
        #
        #     self._debug_collect_state_populations(replicas_to_exchange)
        # self._debug_validate_state_populations(replicas_to_exchange,
        #                                        states_to_exchange,U)
        sampling_time = time.time() - sampling_start_time

    # see children classes for specific implementations
    def update_state_of_replica(self, repl):
        pass

