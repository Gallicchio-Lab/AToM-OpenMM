from __future__ import print_function
"""
Multiprocessing job transport for AsyncRE/OpenMM
"""
import os, re, sys, time, shutil, copy, random, signal
import logging

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from simtk.openmm.app.desmonddmsfile import *
from datetime import datetime
      
class OMMReplica(object):
    #
    # Holds and manages OpenMM system for a replica
    #
    def __init__(self, replica_id, basename):
        self._id = replica_id
        self.basename = basename

        self.pot = None
        self.par = None
        self.is_energy_assigned = False
        self.is_state_assigned = False
        self.cycle = 0
        self.stateid = None
        self.mdsteps = 0
        
        self.open_dms()
        
        self.positions = copy.deepcopy(self.dms.positions)
        self.velocities = copy.deepcopy(self.dms.velocities)        


    def set_state(self, stateid, par):
        self.stateid = int(stateid)
        self.par = par
        self.is_state_assigned = True
        
    def get_state(self):
        return (self.stateid, self.par)

    def get_energy(self):
        return self.pot

    def set_energy(self, pot):
        self.pot = pot
        self.is_energy_assigned = True
        
    def set_posvel(self, positions, velocities):
        self.positions = positions
        self.velocities = velocities

    def open_dms(self):
        input_file  = '%s_0.dms' % self.basename 

        if not os.path.isdir('r%d' % self._id):
            os.mkdir('r%d' % self._id)

        output_file  = 'r%d/%s_ckp.dms' % (self._id,self.basename)
        if not os.path.isfile(output_file):
            shutil.copyfile(input_file, output_file)

        self.dms = DesmondDMSFile([output_file]) 

        self.sql_conn = self.dms._conn[0]
                
        # check for tre_data table in dms file
        tables = self.dms._tables[0]
        conn = self.sql_conn
        if 'tre_data' in tables:
            # read tre_data table
            q = """SELECT epot,temperature,cycle,stateid,mdsteps FROM tre_data WHERE id = 1"""
            ans = conn.execute(q)
            for (epot,temperature,cycle,stateid,mdsteps) in conn.execute(q):
                self.pot = [epot]
                self.par = [temperature]
                self.cycle = cycle
                self.stateid = stateid
                self.mdsteps = mdsteps

    def save_dms(self):
        if self.is_state_assigned and self.is_energy_assigned:
            conn = self.sql_conn
            tables = self.dms._tables[0]
            if 'tre_data' not in tables:
                conn.execute("CREATE TABLE IF NOT EXISTS tre_data (id INTEGER PRIMARY KEY, epot REAL, temperature REAL, cycle INTEGER, stateid INTEGER, mdsteps INTEGER )")
                conn.execute("INSERT INTO tre_data (epot,temperature,cycle,stateid,mdsteps) VALUES (0,0,0,0,0)")
                conn.commit()
                self.dms._tables[0] = self.dms._readSchemas(conn)
            pot_energy =  float(self.pot[0])            
            temperature = float(self.par[0])
            conn.execute("UPDATE tre_data SET epot = %f, temperature = %f, cycle = %d, stateid = %d, mdsteps = %d WHERE id = 1" % (pot_energy, temperature, self.cycle, self.stateid, self.mdsteps))
            conn.commit()
            self.dms.setPositions(self.positions)
            self.dms.setVelocities(self.velocities)

    def set_mdsteps(self, mdsteps):
        self.mdsteps = mdsteps

    def get_mdsteps(self):
        return self.mdsteps

    def set_cycle(self, cycle):
        self.cycle = cycle
        
    def get_cycle(self):
        return self.cycle

    def get_stateid(self):
        return self.stateid

    def _getOpenMMData(self,outfile):
        """
        Reads all of the Openmm simulation data values temperature, energies,
        etc.
        """
        data = []
        f = open(outfile, "r")
        line = f.readline()
        while line:
            datablock = []
            for word in line.split():
                datablock.append(float(word))
            data.append(datablock)
            line = f.readline()
        f.close()
        return data
    
    def set_statepot_from_outputfile(self, replica, cycle):
        outfile = "r%d/%s_%d.out" % (replica, self.basename, cycle)
        data = self._getOpenMMData(outfile)
        # format is temperature, pot_energy
        nr = len(data)
        temperature = float(data[nr-1][0])
        pot_energy = float(data[nr-1][1])
        self.set_state(self.stateid, [temperature])
        self.set_energy([pot_energy])
        
    def set_posvel_from_file(self, replica, cycle):
        tfile = "r%d/%s_%d.dms" % (replica, self.basename, cycle)
        dms = DesmondDMSFile([tfile])
        self.positions = copy.deepcopy(dms.positions)
        self.velocities = copy.deepcopy(dms.velocities)
        dms.close()

    def write_posvel_to_file(self, cyle):
        tfile = "r%d/%s_%d.dms" % (replica, self.basename, cycle)
        dms = DesmondDMSFile([tfile])
        dms.setPositions(self.positions)
        dms.setVelocities(self.velocities)
        dms.close()
