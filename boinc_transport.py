"""
BOINC job transport for AsyncRE
"""
import os
import sys
import time
import re
import pickle
import MySQLdb
import logging
import subprocess

from transport import Transport

class boinc_transport(Transport):
    """
    Class to launch and monitor jobs through a BOINC project
    """
    def __init__(self, jobname, keywords, nreplicas, files_to_stage):
        # jobname: identifies current asyncRE job
        # files_to_stage: permanent files to stage into BOINC download directory (main struct file, datafiles, etc.)
        Transport.__init__(self)
        self.logger = logging.getLogger("async_re.boinc_transport")

        self.jobname = jobname

        # variables required for boinc-based transport
        if keywords.get('BOINC_PROJECTDIR') is None:
            self._exit("BOINC transport requires a BOINC_PROJECTDIR")
        self.project_dir = keywords.get('BOINC_PROJECTDIR')
        if not os.path.isdir(self.project_dir):
            self.logger.critical("Unable to located BOINC project directory: %s",
                                 self.project_dir)
            sys.exit(1)

        #sets up a mapping between replicas running and work unit ids
        self.replica_to_wuid = [ None for k in range(nreplicas)]

        #sets up lookup table for replica done status
        self.replica_status = dict()

        #set connection with mysql boinc server database
        if keywords.get('BOINC_DATABASE') is None:
            self.logger.critical("BOINC transport requires a BOINC_DATABASE")
            sys.exit(1)
        if keywords.get('BOINC_DATABASE_USER') is None:
            self.logger.critical("BOINC transport requires a BOINC_DATABASE_USER")
            sys.exit(1)
        if keywords.get('BOINC_DATABASE_PASSWORD') is None:
            self.logger.critical("BOINC transport requires a BOINC_DATABASE_PASSWORD")
            sys.exit(1)
        self.db_name = keywords.get('BOINC_DATABASE')
        self.db_user = keywords.get('BOINC_DATABASE_USER')
        self.db_pwd = keywords.get('BOINC_DATABASE_PASSWORD')

        # stage files
        if files_to_stage is not None:
            for file in files_to_stage:
                filepath = os.getcwd() + "/" + file
                self._stage_file(filepath,file)

    def restart(self):
        # read from saved file
        self.logger.info("Reading from saved file")
        status_file = "%s_boinc.stat" % self.jobname
        try:
            f = open(status_file,'r')
            self.replica_to_wuid = pickle.load(f)
            for wuid in self.replica_to_wuid:
                if not wuid: continue
                self.replica_status[wuid] = False
            f.close()
        except:
            None

    def save_restart(self):
        #write to saved file
        status_file = "%s_boinc.stat" % self.jobname
        try:
            f = open(status_file,'w')
            pickle.dump(self.replica_to_wuid, f)
            f.close()
        except:
            None

    def _stage_file(self, file_src, file_dest):
        command = "cd %s ; bin/stage_file_v2 %s %s" % (self.project_dir, file_src, file_dest)
        self.logger.info(command)
        res = os.system("cd %s ; bin/stage_file_v2 %s %s" % (self.project_dir, file_src, file_dest))
        if not res == 0:
            self.logger.error("_stage_file(): Error in staging file %s into %s",
                              file_src, file_dest)


    def launchJob(self, replica, job_info):
        """
Enqueues a job based on provided job info.
        """

        if self.replica_to_wuid[replica] != None:
           # a wuid for this replica already exists
           # try to reset it through isDone() ...
           cycle = job_info["cycle"]
           self.isDone(replica,cycle)

        input_file = job_info["input_file"]
        executable = job_info["executable"]
        cycle = job_info["cycle"]
        working_directory = job_info["working_directory"]
        command = "cd %s ; %s %s %s %s %s" % (self.project_dir, executable, working_directory, self.jobname, replica, cycle)

        self.logger.info(command)
        proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        (out, err) = proc.communicate()
        # out should have the workunit name
        try:
            wuid = int(out.split()[0])
        except:
            self.logger.warning("launchjob(): unable to create workunit")
            return None

        old_wuid = self.replica_to_wuid[replica]
        if old_wuid:
            self.logger.info("No longer tracking wuid: %s", old_wuid)
            self.replica_status.pop(old_wuid)

        self.replica_to_wuid[replica] = wuid
        self.replica_status[wuid] = False
        self.logger.info("Now tracking wuid: %s", wuid)

        #write status file
        self.save_restart()

        return 1

    def poll(self):
        self.logger.info("Polling BOINC DB")

        wuids = self.replica_status.keys()
        wuid_strings = map(str, wuids)
        if not wuids:
            self.logger.info("Polling BOINC DB complete! Didn't find any wuids, though")
            return

        self.boinc_db = MySQLdb.connect(user=self.db_user, passwd=self.db_pwd,
                                        db=self.db_name)
        if not self.boinc_db:
            self.logger.critical("poll(): Unable to open connection with mysql db %s", self.db_name)
            sys.exit(1)

        cur = self.boinc_db.cursor()
        sql_command = ("SELECT workunitid, server_state FROM result where "
                       "result.workunitid IN (%s)")
        sql_command = sql_command % ','.join(['%s']*len(wuids))

        cur.execute(sql_command, wuid_strings)

        updated = set()
        for wuid, status in cur.fetchall():
            self.replica_status[wuid] = (status == 5)
            updated.add(wuid)
        cur.close()

        for wuid in set(wuids) - updated:
            self.logger.warning("poll(): cannot locate wu %s in db", str(wuid))
            self.logger.warning("Assume it is not done and hope for the best")

        #test for assimilation
        cur = self.boinc_db.cursor()
        sql_command = ("SELECT id, assimilate_state FROM workunit "
                       "WHERE workunit.id IN (%s)")
        sql_command = sql_command % ','.join(['%s']*len(wuids))
        cur.execute(sql_command, wuid_strings)
        result = cur.fetchall()
        cur.close()

        for wuid, assimilate_state in result:
            if assimilate_state == 2:
                time.sleep(3)
            self.replica_status[wuid] = (assimilate_state == 2)

        self.boinc_db.close()
        self.logger.info("Polling BOINC DB complete!")

    def ProcessJobQueue(self, mintime, maxtime):
        """
with boinc there's no queue to process. Just wait until maxtime
        """
        time.sleep(maxtime)

    def isDone(self,replica,cycle):
        """
          Checks if a replica completed a run.
          Inspects the 'result' mysql boinc server database and considers the workunit complete
          if server_status is 'OVER' (=5).
        """
        wuid = self.replica_to_wuid[replica]

        if wuid == None :
            #unknown work unit, assume it is done and let the file
            #checker discover if it needs to be resubmitted
            return True
        else:
            return self.replica_status[wuid]

if __name__ == "__main__":
    pass
