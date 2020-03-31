"""
SSH job transport for AsyncRE
"""
from __future__ import print_function
from __future__ import division
import os
import re
import sys
import time
import random
import multiprocessing as mp
import logging
import subprocess
import Queue

from ommreplica import OMMReplica
from transport import Transport

class ssh_transport(Transport):
    """
    Class to launch and monitor jobs on a set of nodes via ssh (paramiko)
    """
    def __init__(self, jobname, compute_nodes, openmm_replicas):
        # jobname: identifies current asyncRE job
        # compute_nodes: list of names of nodes in the pool
        # openmm_replicas: list of openmm replica objects
        Transport.__init__(self)
        self.logger = logging.getLogger("async_re.ssh_transport")

        # names of compute nodes (slots)
        self.compute_nodes = compute_nodes
        self.nprocs = len(self.compute_nodes)
        #self.openmm_platform = None

        # node status = None if idle
        # Otherwise a structure containing:
        #    replica number being executed
        #    process id
        #    process name
        #    ...
        self.node_status = [ None for k in range(self.nprocs)]
        #self.openmm_platform = None

        #records replica OpenMM objects
        self.openmm_replicas = openmm_replicas
        nreplicas = len(openmm_replicas)

        # contains the nodeid of the node running a replica
        # None = no information about where the replica is running
        self.replica_to_job = [ None for k in range(nreplicas) ]

        # implements a queue of jobs from which to draw the next job
        # to launch
        self.jobqueue = Queue.Queue()

    def _clear_resource(self, replica):
        # frees up the node running a replica identified by replica id
        job = None
        try:
            job = self.replica_to_job[replica]
        except:
            self.logger.warning("clear_resource(): unknown replica id %d",
                                replica)

        if job == None:
            return None

        try:
            nodeid = job['nodeid']
        except:
            self.logger.warning("clear_resource(): unable to query nodeid")
            return None

        try:
            self.node_status[nodeid] = None
        except:
            self.logger.warning("clear_resource(): unknown nodeid %", nodeid)
            return None

        return nodeid

    def _availableNode(self):
        #returns a node at random among available nodes
        available = [node for node in range(self.nprocs)
                     if self.node_status[node] == None]
        if available == None or len(available) == 0:
            return None
        random.shuffle(available)
        return available[0]

    def _launchCmd(self, command, job):
        if job['username'] =="":
            ssh_command = "ssh %s" % job['nodename']
            scp_remote = "%s:" % job['nodename']
        else:
            ssh_command = "ssh %s@%s" % (job['nodename'],job['username'])
            scp_remote = "%s@%s:" % (job['username'],job['username'])

        if job["remote_working_directory"]:
            mkdir_command = ssh_command + " mkdir -p %s" % job['remote_working_directory']
            self.logger.info(mkdir_command)
            subprocess.call(mkdir_command, shell=True)
            for filename in job["job_input_files"]:
                local_file = job["working_directory"] + "/" + filename
                remote_file = job["remote_working_directory"] + "/" + filename
                send_file_command = "scp " + local_file + " " + scp_remote + remote_file
                self.logger.info(send_file_command)
                subprocess.call(send_file_command, shell=True)


            chmod_command = ssh_command + " chmod -R 777 %s" % job['remote_working_directory']
            self.logger.info(chmod_command)
            subprocess.call(chmod_command, shell=True)

            inputfile = job["remote_working_directory"]+ "/" + job['input_file']
            self.logger.info( "Inputfile: %s Platform: %s Slot number:  %s", inputfile, job['platform'],job['nslots'])
            cfmod_command = ssh_command + " sed -i 's/@platform@/%s/g' %s ;" % (job['platform'],inputfile)
            if job['nslots']:
                cf1mod_command = ssh_command + " sed -i 's/@pn@/%s/g' %s" % (job['nslots'],inputfile)
            else:
                cf1mod_command = ssh_command + " sed -i 's/@pn@//g' %s" % (job['nslots'],inputfile)                
            cfmod_command = cfmod_command + cf1mod_command
            self.logger.info(cfmod_command)
            subprocess.call(cfmod_command, shell=True)

        launch_command = ssh_command + " " + '"%s"' % command
        self.logger.info(launch_command)
        subprocess.call(launch_command, shell=True)

        if job["remote_working_directory"]:
            for filename in job["job_output_files"]:
                local_file = job["working_directory"] + "/" + filename
                remote_file = job["remote_working_directory"] + "/" + filename
                try:
                    receive_file_command = "scp " + scp_remote + remote_file + " " + local_file
                    subprocess.call(receive_file_command, shell=True)
                except:
                    self.logger.info("Warning: unable to copy back file %s" % local_file)

            rmdir_command = ssh_command + " rm -rf %s" % job['remote_working_directory']
            subprocess.call(rmdir_command, shell=True)


    def launchJob(self, replica, job_info):
        """
        Enqueues a job based on provided job info.
        """
        input_file = job_info["input_file"]
        output_file = job_info["output_file"]
        error_file = job_info["error_file"]
        executable = job_info["executable"]

        command = "%s %s > %s 2> %s " % ( executable, input_file, output_file, error_file)

        job = job_info
        job['replica'] = replica
        job['command'] = command
        job['process_handle'] = None
        job['start_time'] = 0

        self.replica_to_job[replica] = job

        self.jobqueue.put(replica)

        return self.jobqueue.qsize()

    #intel coprocessor setup
    #edit on 10.20.15
    def ModifyCommand(self,job, command):
        nodename = job['nodename']
        nodeN = job['nthreads']
        slotN = job['nslots']
        architecture = job['architecture']

        #add command to go to remote working directory
        cd_to_command = "cd %s ; " % job["remote_working_directory"]

        #mic_pattern = re.compile("mic0" or "mic1")

        #if re.search(mic_pattern, nodename):
        #    offset = slotN * (nodeN/4)
        #    add_to_command = "export KMP_PLACE_THREADS=6C,4T,%dO ; " % offset
        #else:
        #    add_to_command = "export OMP_NUM_THREADS=%d;"% nodeN
        #new_command = add_to_command + cd_to_command + command

        #set up the python env environment 10.21.15
        new_command = cd_to_command + command
        self.logger.info(new_command) #can print new_command here to check the command
        return new_command
    #edit end on 10.20.15

    def DrainJobQueue(self):
        # not needed for ssh transport (?)
        pass

    def ProcessJobQueue(self, mintime, maxtime):
        """
        Launches jobs waiting in the queue.
        It will scan free nodes and job queue up to maxtime.
        If the queue becomes empty, it will still block until maxtime is elapsed.
        """
        njobs_launched = 0
        usetime = 0
        nreplicas = len(self.replica_to_job)

        while usetime < maxtime:

            # find an available node
            node = self._availableNode()

            while (not self.jobqueue.empty()) and (not node == None):

                # grabs job on top of the queue
                replica = self.jobqueue.get()
                job = self.replica_to_job[replica]

                # assign job to available node
                job['nodeid'] = node
                job['nodename'] = self.compute_nodes[node]["node_name"]
                job['nthreads'] = int(self.compute_nodes[node]["threads_number"]) 
                #job['nslots']   = int(self.compute_nodes[node]["slot_number"])
                job['nslots']   = self.compute_nodes[node]["slot_number"]
                job['username'] = self.compute_nodes[node]["user_name"]
                job['openmm_replica'] = self.openmm_replicas[replica]
                #added on 10.22.15 save the arch information for platform, save the slot information for multiple GPU
                job['architecture'] = self.compute_nodes[node]["arch"]
                opencl = 'OpenCL'
                if opencl in job['architecture']:
                    job['platform']='OpenCL'
                else:
                    job['platform']='Reference'
                #added end on 10.22.15

                # get the shell command
                command = job['command']
                #retrieve remote working directory of node
                job["remote_working_directory"] = self.compute_nodes[node]["tmp_folder"] + "/" + job["remote_replica_dir"]

                command=self.ModifyCommand(job,command)

                if job["remote_working_directory"] and job['job_input_files']:
                    for filename in job['job_input_files']:
                        local_file = job["working_directory"] + "/" + filename
                        remote_file = job["remote_working_directory"] + "/" + filename
                        #self.logger.info("%s %s", local_file, remote_file) #can print out here to verify files

                #if self.compute_nodes[node]["arch"]:
                #    architecture = self.compute_nodes[node]["arch"] 
                #self.platforms = architecture
                #else:
                #    architecture = "Reference"
                #self.openmm_platform = architecture
                #self.openmm_slot = job['nslots']

                #add the architecture information and the slot information into the input file #10.22.15
                #inpfile = "r%d/%s_%d.py" % (replica, basename, cycle)

                #end on 10.22.15


                #edit on 10.20.15
                #exec_directory = job["exec_directory"]
                #lib_directory = exec_directory + "/lib/" + architecture
                #bin_directory  = exec_directory + "/bin/" + architecture

                #job["exec_files"] = []
                #for filename in os.listdir(lib_directory):
                #    job["exec_files"].append(lib_directory + "/" + filename)
                #for filename in os.listdir(bin_directory):
                #    job["exec_files"].append(bin_directory + "/" + filename)
                #edit end on 10.20.15
                # launches job
                processid = mp.Process(target=self._launchCmd, args=(command, job))
                processid.start()

                job['process_handle'] = processid

                job['start_time'] = time.time()

                # connects node to replica
                self.replica_to_job[replica] = job
                self.node_status[node] = replica

                # updates number of jobs launched
                njobs_launched += 1
                node = self._availableNode()

            # waits mintime second and rescans job queue
            time.sleep(mintime)

            # updates set of free nodes by checking for replicas that have exited
            for repl in range(nreplicas):
                self.isDone(repl,0)

            usetime += mintime

        return njobs_launched

    def cancel(self, replica):
        """
        Kills a running replica
        """
        job = self.replica_to_job[replica]
        if job == None:
            return

        process = job['process_handle']
        if process == None:
            return
        process.terminate()


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
            process = job['process_handle']
            if process == None:
                done = False
            else:
                done = not process.is_alive()
            if done:
                # update openmm replica
                thiscycle = int(job["cycle"])
                thisreplica = int(job["replica"])
                ommreplica = job['openmm_replica']
                ommreplica.set_cycle(thiscycle+1)
                #attempt retrieve data from output files
                try:
                    ommreplica.set_statepot_from_outputfile(thisreplica,thiscycle)
                    ommreplica.set_posvel_from_file(thisreplica,thiscycle)
                except:
                    self.logger.warning("Unable to retrieve data for replica %d, cycle %d" % (replica,thiscycle))
                # disconnects replica from job and node
                self._clear_resource(replica)
                self.replica_to_job[replica] = None
            elif process:
                time_interval = time.time() - job['start_time']
                if time.time() - job['start_time'] > 18000:
                    self.logger.info("time interval is %f for replica %d", time_interval, replica)
                    self.logger.info("18000 seconds time limit is exceeded for replica %d", replica)
                    self.cancel(replica)
                    self._clear_resource(replica)
                    self.replica_to_job[replica] = None
                    done = True
            return done
