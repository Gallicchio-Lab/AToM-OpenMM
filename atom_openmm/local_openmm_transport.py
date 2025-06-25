from __future__ import print_function
from __future__ import division
"""
Multiprocessing job transport for AsyncRE/OpenMM
"""
import os, re, sys, time, shutil, copy, random, signal
import multiprocessing as mp
#from multiprocessing import Process, Queue, Event
import logging

import openmm as mm
from openmm.app import *
from openmm import *
from openmm.unit import *
from datetime import datetime

from atom_openmm.ommreplica import *
from atom_openmm.ommsystem import *
from atom_openmm.ommworker import *

from contextlib import contextmanager

from atom_openmm.transport import Transport

class LocalOpenMMTransport(Transport):
    """
    Class to launch and monitor jobs on a set of local GPUs
    """
    def __init__(self, jobname, openmm_workers, openmm_replicas):
        # jobname: identifies current asyncRE job
        Transport.__init__(self)
        self.logger = logging.getLogger("async_re.local_openmm_transport")

        # openmm contexts
        self.openmm_workers = openmm_workers
        self.nprocs = len(self.openmm_workers)

        # record replica OpenMM objects
        self.openmm_replicas = openmm_replicas

        # device status = None if idle
        # Otherwise a structure containing:
        #    replica id being executed
        #    ...
        self.node_status = [ None for k in range(self.nprocs)]

        # contains information about the device etc. running a replica
        # None = no information about where the replica is running
        self.replica_to_job = [ None for k in range(len(openmm_replicas)) ]

        # implements a queue of jobs from which to draw the next job
        # to launch
        ctx = mp.get_context('spawn')
        self.jobqueue = ctx.Queue()

        self.ncrashes = [ 0 for k in range(self.nprocs)]
        self.disabled = [ False for k in range(self.nprocs)]
        self.maxcrashes = 4
        self.done_replicas = []

    def _clear_resource(self, replica):
        # frees up the node running a replica identified by replica id
        job = {}
        try:
            job = self.replica_to_job[replica]
            if job is None:
                return None
        except:
            self.logger.warning("clear_resource(): unknown replica id %d",
                                replica)

        if 'nodeid' not in job:
            return None
        else:
            nodeid = job['nodeid']

        try:
            if self.node_status[nodeid] is not None and self.node_status[nodeid] >= 0: #-1 signals a crashed node that should be left alone
                self.node_status[nodeid] = None
        except:
            self.logger.warning("clear_resource(): unable to query nodeid %d", nodeid)
            return None

        return nodeid

    def numNodesAlive(self):
        alive = [node for node in range(self.nprocs)
                     if self.node_status[node] is None or self.node_status[node] >= 0 ]
        return len(alive)

    def fixnodes(self):
        for nodeid in range(self.nprocs):
            if self.node_status[nodeid] is not None and self.node_status[nodeid] < 0 and not self.disabled[nodeid]:
                if self.ncrashes[nodeid] <= self.maxcrashes:
                    self.ncrashes[nodeid] += 1
                    self.logger.warning("fixnodes(): attempting to restart nodeid %d", nodeid)
                    res = self.openmm_workers[nodeid].start_worker()
                    if res is not None:
                        self.node_status[nodeid] = None
                else:
                    self.logger.warning("fixnodes(): node %d has crashed too many times; it will not be restarted.", nodeid)
                    self.disabled[nodeid] = True

    def _availableNode(self):
        #returns a node at random among available nodes
        available = [node for node in range(self.nprocs)
                     if self.node_status[node] is None]

        if available == None or len(available) == 0:
            return None
        random.shuffle(available)

        return available[0]

    def launchJob(self, replica, job_info):
        #Enqueues a replica for running based on provided job info.
        job = job_info
        job['replica'] = replica
        job['start_time'] = 0
        self.replica_to_job[replica] = job
        self.jobqueue.put(replica)
        return self.jobqueue.qsize()

    def LaunchReplica(self, worker, replica, cycle, nsteps,
                      nheating = 0, ncooling = 0, hightemp = 0.0):
        #if replica.contextchkpt != None:
        #    worker.set_chkpt(replica.contextchkpt)
        (stateid, par) = replica.get_state()
        worker.set_posvel(replica.positions, replica.velocities)
        worker.set_state(par)
        worker.run(nsteps, nheating, ncooling, hightemp)

    def ProcessJobQueue(self, mintime, maxtime):
        #Launches jobs waiting in the queue.
        #It will scan free devices and job queue up to maxtime.
        #If the queue becomes empty, it will still block until maxtime is elapsed.
        njobs_launched = 0
        nreplicas = len(self.replica_to_job)

        when_started = time.time()
        while time.time() < when_started + maxtime:
            # find an available node
            node = self._availableNode()
            while (not self.jobqueue.empty()) and (not node == None):

                # grabs job on top of the queue
                replica = self.jobqueue.get()
                job = self.replica_to_job[replica]

                # assign job to available node
                job['nodeid'] = node
                job['openmm_replica'] = self.openmm_replicas[replica]
                job['openmm_worker'] = self.openmm_workers[node]
                job['start_time'] = time.time()

                # connects node to replica
                self.replica_to_job[replica] = job
                self.node_status[node] = replica

                if 'nheating' in job:
                    nheating = job['nheating']
                    ncooling = job['ncooling']
                    hightemp = job['hightemp']
                else:
                    nheating = 0
                    ncooling = 0
                    hightemp = 0.0

                self.LaunchReplica(job['openmm_worker'], job['openmm_replica'], job['cycle'],
                                   job['nsteps'], nheating, ncooling, hightemp)

                # updates number of jobs launched
                njobs_launched += 1
                if not hasattr(job['openmm_worker'], '_runningSignal'):
                    # Synchronous single worker. Set node_status to None to reuse for next replica.
                    self.node_status[node] = None
                node = self._availableNode()

            self._finalize_replicas()

            # waits mintime second and rescans job queue
            time.sleep(mintime)

            # updates set of free nodes by checking for replicas that have exited
            for repl in range(nreplicas):
                self.isDone(repl,0)

            #restarts crashed nodes if any
            self.fixnodes()

        self._finalize_replicas()

        return njobs_launched

    def DrainJobQueue(self):
        #clear the job queue
        while not self.jobqueue.empty():
            # grabs job on top of the queue
            replica = self.jobqueue.get()
            self._clear_resource(replica)
            self.replica_to_job[replica] = None

    def _update_replica(self, job):
        import numpy as np

        #update replica cycle, mdsteps, write out, etc. from worker
        ommreplica = job['openmm_replica']
        if job['openmm_worker'].has_crashed(): #refuses to update replica from a crashed worker
            return None
        (pos,vel) = job['openmm_worker'].get_posvel()
        #chkpt = job['openmm_worker'].get_chkpt()
        pot = job['openmm_worker'].get_energy()
        if pos is None or vel is None or pot is None:
            return None
        for value in pot.values():
            if math.isnan(value._value):
                return None
        if np.any(np.isnan(pos)):
            return None
        if np.any(np.isnan(vel)):
            return None
        cycle = ommreplica.get_cycle() + 1
        ommreplica.set_cycle(cycle)
        mdsteps = ommreplica.get_mdsteps() + job['nsteps']
        ommreplica.set_mdsteps(mdsteps)
        #update positions and velocities of openmm replica
        ommreplica.set_posvel(pos,vel)
        #ommreplica.set_chkpt(chkpt)
        
        #TODO: should also update boxsize
        #update energies of openmm replica
        ommreplica.set_energy(pot)
        #output data and trajectory file update 
        if mdsteps % job['nprnt'] == 0:
            ommreplica.save_out()
        #queue replica to perform time-consuming tasks such as trajectory writing to
        #after new replicas are launched
        self.done_replicas.append((ommreplica, mdsteps, job['ntrj']))
        return 0

    def _finalize_replicas(self):
        for donereplica in self.done_replicas:
            ommreplica = donereplica[0]
            mdsteps = donereplica[1]
            ntrj = donereplica[2]
            if mdsteps % ntrj == 0:
                ommreplica.save_xtc()
        self.done_replicas = []

    def isDone(self,replica,cycle):
        """
        Checks if a replica completed a run.

        If a replica is done it clears the corresponding node.
        Note that cycle is ignored by job transport. It is assumed that it is
        the latest cycle.  it's kept for argument compatibility with
        hasCompleted() elsewhere.
        """
        job = self.replica_to_job[replica]
        if job == None:
            # if job has been removed we assume that the replica is done
            return True
        else:
            try:
                openmm_worker = job['openmm_worker']
            except:
                #job is in the queue but not yet launched
                return False

            if openmm_worker.has_crashed():
                self.logger.warning("isDone(): replica %d has crashed", replica)
                openmm_worker.finish(wait = False)
                self.node_status[job['nodeid']] = -1 #signals dead context
                self._clear_resource(replica)
                self.replica_to_job[replica] = None
                return True

            if not openmm_worker.is_started():
                done = False
            else:
                done = openmm_worker.is_running() and openmm_worker.is_done()

            if done:
                #update replica info
                if hasattr(openmm_worker, '_runningSignal'):
                    openmm_worker._runningSignal.clear()
                retcode = self._update_replica(job)
                if retcode is None:
                    self.logger.warning("isDone(): replica %d has completed with errors", replica)
                    self.node_status[job['nodeid']] = -1 #signals dead context
                # disconnects replica from job and node
                self._clear_resource(replica)
                #flag replica as not linked to a job
                self.replica_to_job[replica] = None

            return done
