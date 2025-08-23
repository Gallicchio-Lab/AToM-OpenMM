from __future__ import print_function
"""
Multiprocessing job transport for AsyncRE/OpenMM
"""
import os
import copy
from openmm.unit import kelvin, kilojoules_per_mole, kilocalories_per_mole
from openmm.app import XTCFile
from openmm import unit

class OMMReplica(object):
    #
    # Holds and manages OpenMM state for a replica
    #
    def __init__(self, replica_id, basename, worker, logger, keywords):
        self._id = replica_id
        self.basename = basename
        self.worker = worker
        self.context = worker.context
        self.ommsystem = worker.ommsystem
        self.logger = logger
        self.keywords = keywords
        self.pot = None
        self.par = None
        self.cycle = 1
        self.stateid = None
        self.mdsteps = 0
        self.outfile = None
        self.safeckpt_file = "ckpt_is_valid"

        state = self.context.getState(getPositions=True, getVelocities=True)
        self.positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
        self.velocities = state.getVelocities(asNumpy=True).value_in_unit(unit.nanometers/unit.picoseconds)
        self.contextchkpt = None #holds a complete state of the Context
        
        if not os.path.isdir('r%d' % self._id):
            os.mkdir('r%d' % self._id)
            
        self.open_out()
        #may override stateid, positions, etc.
        self.load_checkpoint()
        self.open_xtc()

    def set_state(self, stateid, par):
        self.stateid = int(stateid)
        self.par = copy.deepcopy(par)
        self.is_state_assigned = True
        self.update_context_from_state()
        
    def get_state(self):
        return (self.stateid, self.par)

    def get_energy(self):
        return self.pot

    def set_energy(self, pot):
        self.pot = copy.deepcopy(pot)
        
    def set_posvel(self, positions, velocities):
        self.positions[:] = positions
        self.velocities[:] = velocities

    def set_chkpt(self, chkpt):
        self.contextchkpt = chkpt

    def open_out(self):
        outfilename =  'r%d/%s.out' % (self._id,self.basename)
        self.outfile = open(outfilename, 'a+')
        if self.outfile is None:
            self.logger.warning("unable to open outfile %s" % outfilename)

    def load_checkpoint(self):
        override_safecheckpoint = self.keywords.get('OVERRIDE_SAFECHECKPOINT')
        ckptname = self.keywords.get('CHECKPOINT_FILE', f'{self.basename}_ckpt.xml')
        ckptfile = 'r%d/%s' % (self._id,ckptname)
        if os.path.isfile(ckptfile):
            if os.path.isfile(self.safeckpt_file) or (override_safecheckpoint is not None):
                self.logger.info("Loading checkpointfile %s" % ckptfile) 
                self.worker.simulation.loadState(ckptfile)
                self.update_state_from_context()
            else:
                self.logger.error("The simulation has been interrupted while writing the checkpoint files. The checkpoint files might be corrupted. Remove the replica directories and restart from scratch. Alternatively, force the loading of the checkpoint files by setting OVERRIDE_SAFECHECKPOINT in the control file.")
                raise ValueError('Bad checkpoints')

    def save_checkpoint(self):
        if not os.path.isfile(self.safeckpt_file):#refuse to write checkpoint in unsafe mode
            ckptname = self.keywords.get('CHECKPOINT_FILE', f'{self.basename}_ckpt.xml')
            ckptfile = 'r%d/%s' % (self._id,ckptname)
            self.update_context_from_state()
            self.worker.simulation.saveState(ckptfile)
        else:
           self.logger.warning("Refused attempt to save checkpoint file %s in unsafe mode. Remove file %s prior to writing checkpoints and restore it when done." % (ckptfile, self.safeckpt_file) )

    def open_xtc(self):
        xtcfilename =  'r%d/%s.xtc' % (self._id,self.basename)
        interval = int(self.keywords.get('TRJ_FREQUENCY'))
        append = os.path.isfile(xtcfilename)
        self.xtcfile = xtcfilename
        self.xtc = XTCFile(self.xtcfile, self.worker.topology, self.ommsystem.MDstepsize, interval=interval, append=append)


    def save_xtc(self):
        #TODO
        #boxsize options works only for NVT because the boxsize of the service worker
        #is not updated from the compute worker
        boxsize = self.worker.simulation.context.getState().getPeriodicBoxVectors()
        self.xtc.writeModel(self.positions, periodicBoxVectors=boxsize)

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


class OMMReplicaTRE(OMMReplica):
    def save_out(self):
        if self.pot is not None and self.par is not None:
            pot_energy = self.pot['potential_energy']
            temperature = self.par['temperature']
            if self.outfile is not None:
                self.outfile.write("%d %f %f\n" % (self.stateid, temperature, pot_energy))

    def update_state_from_context(self):
        self.cycle = int(self.context.getParameter(self.ommsystem.parameter['cycle']))
        self.stateid = int(self.context.getParameter(self.ommsystem.parameter['stateid']))
        self.mdsteps = int(self.context.getParameter(self.ommsystem.parameter['mdsteps']))
        if self.par is None:
            self.par = {}
        self.par['temperature'] = self.context.getParameter(self.ommsystem.parameter['temperature'])*kelvin
        if self.pot is None:
            self.pot = {}
        self.pot['potential_energy'] = self.context.getParameter(self.ommsystem.parameter['potential_energy'])*kilojoules_per_mole
        state = self.context.getState(getPositions=True, getVelocities=True)
        self.positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
        self.velocities = state.getVelocities(asNumpy=True).value_in_unit(unit.nanometers/unit.picoseconds)

    def update_context_from_state(self):
        self.context.setParameter(self.ommsystem.parameter['cycle'], self.cycle)
        self.context.setParameter(self.ommsystem.parameter['stateid'], self.stateid)
        self.context.setParameter(self.ommsystem.parameter['mdsteps'], self.mdsteps)
        if self.par is None:
            self.context.setParameter(self.ommsystem.parameter['temperature'], self.par['temperature']/kelvin)
        if self.pot is None:
            self.context.setParameter(self.ommsystem.parameter['potential_energy'], self.pot['potential_energy']/kilojoules_per_mole)

class OMMReplicaATM(OMMReplica):
    def get_output_data(self):
        # return the output data for current state as a tuple
        if self.pot is not None and self.par is not None:
            pot_energy = self.pot["potential_energy"]
            pert_energy = self.pot["perturbation_energy"]
            bias_energy = self.pot["bias_energy"]
            temperature = self.par["temperature"]
            lmbd1 = self.par["lambda1"]
            lmbd2 = self.par["lambda2"]
            alpha = self.par["alpha"]
            uh = self.par["uh"]
            w0 = self.par["w0"]
            direction = self.par["atmdirection"]

            replica_state_data = (
                self.stateid,
                temperature / kelvin,
                direction,
                lmbd1,
                lmbd2,
                alpha * kilocalories_per_mole,
                uh / kilocalories_per_mole,
                w0 / kilocalories_per_mole,
                pot_energy / kilocalories_per_mole,
                pert_energy / kilocalories_per_mole,
                bias_energy / kilocalories_per_mole,
            )

            return replica_state_data
        else:
            self.logger.error("unable to save output")
            raise ValueError
    
    def save_out(self, data=None):
        if data is None:
            try:
                data = [self.get_output_data()]
            except ValueError:
                self.logger.error("unable to save output")

        if self.outfile is not None:
            for line in data:
                self.outfile.write("%d %f %f %f %f %f %f %f %f %f %f\n" % line)
            self.outfile.flush()
        else:
            self.logger.warning("unable to save output")
        
    def update_state_from_context(self):
        self.cycle = int(self.context.getParameter(self.ommsystem.parameter['cycle']))
        self.stateid = int(self.context.getParameter(self.ommsystem.parameter['stateid']))
        self.mdsteps = int(self.context.getParameter(self.ommsystem.parameter['mdsteps']))
        if self.par is None:
            self.par = {}
        self.par['temperature'] = self.context.getParameter(self.ommsystem.parameter['temperature'])*kelvin
        self.par['lambda1'] = self.context.getParameter(self.ommsystem.atmforce.Lambda1())
        self.par['lambda2'] = self.context.getParameter(self.ommsystem.atmforce.Lambda2())
        self.par['alpha'] = self.context.getParameter(self.ommsystem.atmforce.Alpha())/kilojoules_per_mole
        self.par['uh'] = self.context.getParameter(self.ommsystem.atmforce.Uh())*kilojoules_per_mole
        self.par['w0'] = self.context.getParameter(self.ommsystem.atmforce.W0())*kilojoules_per_mole
        self.par['atmdirection'] = self.context.getParameter(self.ommsystem.atmforce.Direction())
        self.par['atmintermediate'] = self.context.getParameter(self.ommsystem.parameter['atmintermediate'])
        self.par[self.ommsystem.atmforce.Umax()] = self.context.getParameter(self.ommsystem.atmforce.Umax())*kilojoules_per_mole
        self.par[self.ommsystem.atmforce.Ubcore()] = self.context.getParameter(self.ommsystem.atmforce.Ubcore())*kilojoules_per_mole
        self.par[self.ommsystem.atmforce.Acore()] = self.context.getParameter(self.ommsystem.atmforce.Acore())
        if self.pot is None:
            self.pot = {}
        self.pot['potential_energy'] = self.context.getParameter(self.ommsystem.parameter['potential_energy'])*kilojoules_per_mole
        self.pot['perturbation_energy'] = self.context.getParameter(self.ommsystem.parameter['perturbation_energy'])*kilojoules_per_mole
        self.pot['bias_energy'] = self.context.getParameter(self.ommsystem.parameter['bias_energy'])*kilojoules_per_mole
        state = self.context.getState(getPositions=True, getVelocities=True)
        self.positions = state.getPositions(asNumpy=True).value_in_unit(unit.nanometers)
        self.velocities = state.getVelocities(asNumpy=True).value_in_unit(unit.nanometers/unit.picoseconds)

    def update_context_from_state(self):
        self.context.setParameter(self.ommsystem.parameter['cycle'], self.cycle)
        self.context.setParameter(self.ommsystem.parameter['stateid'], self.stateid)
        self.context.setParameter(self.ommsystem.parameter['mdsteps'], self.mdsteps)
        if self.par is not None:
            self.context.setParameter(self.ommsystem.parameter['temperature'], self.par['temperature']/kelvin)
            self.context.setParameter(self.ommsystem.atmforce.Lambda1(), self.par['lambda1'])
            self.context.setParameter(self.ommsystem.atmforce.Lambda2(), self.par['lambda2'])
            self.context.setParameter(self.ommsystem.atmforce.Alpha(), self.par['alpha']*kilojoules_per_mole)
            self.context.setParameter(self.ommsystem.atmforce.Uh(), self.par['uh']/kilojoules_per_mole)
            self.context.setParameter(self.ommsystem.atmforce.W0(), self.par['w0']/kilojoules_per_mole)
            self.context.setParameter(self.ommsystem.atmforce.Direction(), self.par['atmdirection'])
            self.context.setParameter(self.ommsystem.parameter['atmintermediate'], self.par['atmintermediate'])
            self.context.setParameter(self.ommsystem.atmforce.Umax(), self.par[self.ommsystem.atmforce.Umax()]/kilojoules_per_mole)
            self.context.setParameter(self.ommsystem.atmforce.Ubcore(), self.par[self.ommsystem.atmforce.Ubcore()]/kilojoules_per_mole)
            self.context.setParameter(self.ommsystem.atmforce.Acore(), self.par[self.ommsystem.atmforce.Acore()])
        if self.pot is not None:
            self.context.setParameter(self.ommsystem.parameter['potential_energy'], self.pot['potential_energy']/kilojoules_per_mole)
            self.context.setParameter(self.ommsystem.parameter['perturbation_energy'], self.pot['perturbation_energy']/kilojoules_per_mole)
            self.context.setParameter(self.ommsystem.parameter['bias_energy'], self.pot['bias_energy']/kilojoules_per_mole)
        self.context.setPositions(self.positions)
        self.context.setVelocities(self.velocities)
