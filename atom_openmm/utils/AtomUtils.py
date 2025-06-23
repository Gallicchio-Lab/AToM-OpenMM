import string
import random
import numpy as np
from contextlib import contextmanager
from pathlib import Path
import os
from openmm import version as ommversion
from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *


@contextmanager
def set_directory(path: Path):
    """Sets the cwd within the context

    Args:
        path (Path): The path to the cwd

    Yields:
        None
    """

    origin = Path().absolute()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)


class AtomUtils(object):
    """
    ATM Force python utilities module
    """
    
    def __init__(self, system, fix_zero_LJparams = True):
        """
        Creates an ATMMetaForceUtils object.

        Parameters
        ----------
        system : System
            The OpenMM System to operate upon
        fix_zero_LJparams: bool
            Whether to fix atoms with zero LJ parameters
        
        Returns
        -------
        An ATMMetaForceUtils object

        """

        self.system = system

        self.CMCMDistForce = None
        self.CMAngleThetaForce = None
        self.CMAnglePhiForce = None
        self.CMAnglePsiForce = None
        self.TorsionalRestraintForce = None

        self.major_ommversion = int(ommversion.version.split(".")[0])
        self.minor_ommversion = int(ommversion.version.split(".")[1])

        if fix_zero_LJparams:
            import re
            nbpattern = re.compile(".*Nonbonded.*")
            for force in self.system.getForces():
                if nbpattern.match(str(type(force))):
                    self.fixZeroLJParams(force)

    def addRestraintForce(self,
                          lig_cm_particles = None, rcpt_cm_particles = None,
                          kfcm = 0.0 * kilocalorie_per_mole/angstrom**2,
                          tolcm = 0.0 * angstrom,
                          offset = [0., 0., 0.] * angstrom):
        """
        Deprecated function. Calls addVsiteRestraintForceCMCM()
        """ 
        print("warning: AddRestraintForce() is deprecated. Use addVsiteRestraintForceCMCM()")
        res = self.addVsiteRestraintForceCMCM(lig_cm_particles, rcpt_cm_particles,
                                              kfcm, tolcm, offset)
        return res

    #
    # Adds restrains to keep the the CM of the ligand atoms around the CM of the receptor atoms
    # plus an offset, within a tolerance, with a flat-bottom harmonic potential
    def addVsiteRestraintForceCMCM(self,
                          lig_cm_particles = None, rcpt_cm_particles = None,
                          kfcm = 0.0 * kilocalorie_per_mole/angstrom**2,
                          tolcm = 0.0 * angstrom,
                          offset = [0., 0., 0.] * angstrom):
        """
        Adds a flat-bottom quadratic potential restraint to keep the the CM of the ligand atoms 
        near the CM of the receptor atoms. This potential is used to define the spatial extent of
        the binding site.
        
        Parameters
        ----------
        lig_cm_particles: list of ints
            List of the indexes of the atoms that define the CM of the ligand molecule
        rcpt_cm_particles : list of ints
            List of the indexes of the atoms that define the CM of the receptor molecule
        kfcm : Quantity in energy/distance^2 units
            Force constant of the quadratic potential
        tolcm : Quantity in distance units
            Tolerance of the flat-bottom restraint
        offset : List of 3 Quantity's in distance units
            Center of the restraint. The potential is applied to the CM-CM distance after subtracting offset.
            It is used to keep the CM-CM distance near offset.

        Returns
        -------
        The CM-CM restraint Force object added to the OpenMM system

        """
        assert len(lig_cm_particles) > 0
        assert len(rcpt_cm_particles) > 0

        expr = "(kfcm/2)*step(d12-tolcm)*(d12-tolcm)^2 "
        expr += " ; d12 = sqrt((x1 - offx - x2)^2 + (y1 - offy - y2)^2 + (z1 - offz - z2)^2 ) ; "
        force = None
        if self.CMCMDistForce is None:
           self.CMCMDistForce  =  mm.CustomCentroidBondForce(2,expr)
           force = self.CMCMDistForce
           force.addPerBondParameter("kfcm")
           force.addPerBondParameter("tolcm")
           force.addPerBondParameter("offx")
           force.addPerBondParameter("offy")
           force.addPerBondParameter("offz")
           numgroups = 0
           self.system.addForce(force)
        else:
            force = self.CMCMDistForce
            numgroups = force.getNumGroups()
        force.addGroup(lig_cm_particles) #g1 CM of lig
        force.addGroup(rcpt_cm_particles) #g2 CM of rcpt
        kfc = kfcm / (kilojoule_per_mole/nanometer**2)
        tolc = tolcm / nanometer
        offv = offset / nanometer
        offx = offv[0]
        offy = offv[1]
        offz = offv[2]
        groups = [numgroups + 0, numgroups + 1]
        params = [kfc, tolc, offx, offy, offz]
        force.addBond(groups, np.array(params, dtype=np.double))
        return force

    def _roundExpression(self, n, x):
        return "{n} = select(step({x} - floor({x})  - 0.5), ceil({x}), floor({x}))".format(n = n, x = x)

    def _wrapExpression(self, w, x, r):
        letters = string.ascii_lowercase
        n = "n" + ''.join(random.choice(letters) for i in range(4))
        expr = "{w} = {x} - {r}*{n} ; ".format(w = w, x = x, r = r, n = n) + " ; "
        expr += self._roundExpression("{n}".format(n = n), "({x}/{r})".format(x = x, r = r))
        return expr

    def _dotExpression(self, dot, x1, y1, z1, x2, y2, z2):
        return "{dot} = {x1}*{x2} + {y1}*{y2} + {z1}*{z2}".format(dot = dot, x1 = x1, y1 = y1, z1 = z1,  x2 = x2, y2 = y2, z2 = z2) 

    def _unitvExpression(self, xu, yu, zu, x, y, z):
        letters = string.ascii_lowercase
        d = "d" + ''.join(random.choice(letters) for i in range(4))
        return "{xu} = {x}/{d} ; {yu} = {y}/{d} ; {zu} = {z}/{d} ; {d} = sqrt({x}*{x}+{y}*{y}+{z}*{z})".format(
            d = d, xu = xu, yu = yu, zu = zu, x = x, y = y, z = z)

    def _diffvExpression(self, dx, dy, dz, x2, y2, z2, x1, y1, z1):
        return "{dx} = {x2} - {x1} ; {dy} = {y2} - {y1} ; {dz} = {z2} - {z1} ".format(
            dx = dx, dy = dy, dz = dz,
            x2 = x2, y2 = y2, z2 = z2,
            x1 = x1, y1 = y1, z1 = z1)

    def _crossExpression(self, xc, yc, zc, x1, y1, z1, x2, y2, z2):
        return ( "{xc} =  {y1}*{z2} - {y2}*{z1} ; "
                 "{yc} = -{x1}*{z2} + {x2}*{z1} ; "
                 "{zc} =  {x1}*{y2} - {x2}*{y1}   " ).format(xc = xc, yc = yc, zc = zc, x1 = x1, y1 = y1, z1 = z1, x2 = x2, y2 = y2, z2 = z2)

    def _cosangleExpression(self, costh, x1, y1, z1, x2, y2, z2, x3, y3, z3):
        if self.major_ommversion >= 7 and self.minor_ommversion >= 6:
            return costh + "= cos(pointangle(" + x1 + "," + y1 + "," + z1 + "," + \
                                                 x2 + "," + y2 + "," + z2 + "," + \
                                                 x3 + "," + y3 + "," + z3 + "))"
        else:
            letters = string.ascii_lowercase
            v1  = "v1"  + ''.join(random.choice(letters) for i in range(4))
            v1u = "v1u" + ''.join(random.choice(letters) for i in range(4))
            v2  = "v2"  + ''.join(random.choice(letters) for i in range(4))
            v2u = "v2u" + ''.join(random.choice(letters) for i in range(4))
            expr  = self._dotExpression(costh, "x"+v1u,"y"+v1u,"z"+v1u,"x"+v2u,"y"+v2u,"z"+v2u) + " ; "
            expr += self._unitvExpression("x"+v1u,"y"+v1u,"z"+v1u, "x"+v1,"y"+v1,"z"+v1) + " ; "
            expr += self._unitvExpression("x"+v2u,"y"+v2u,"z"+v2u, "x"+v2,"y"+v2,"z"+v2) + " ; "
            expr += self._diffvExpression("x"+v1, "y"+v1, "z"+v1, x1, y1, z1, x2, y2, z2) + " ; "
            expr += self._diffvExpression("x"+v2, "y"+v2, "z"+v2, x3, y3, z3, x2, y2, z2)
            return expr

    def _dihedralExpression(self, phi, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4):
        if self.major_ommversion >= 7 and self.minor_ommversion >= 6:
            return phi + "= pointdihedral(" + x1 + "," + y1 + "," + z1 + "," + \
                                              x2 + "," + y2 + "," + z2 + "," + \
                                              x3 + "," + y3 + "," + z3 + "," + \
                                              x4 + "," + y4 + "," + z4 + ")"
        else:
            letters = string.ascii_lowercase
            v1   = "v1"   + ''.join(random.choice(letters) for i in range(4))
            v2   = "v2"   + ''.join(random.choice(letters) for i in range(4))
            v2u  = "v2u"  + ''.join(random.choice(letters) for i in range(4))
            v3   = "v3"   + ''.join(random.choice(letters) for i in range(4))
            n1   = "n1"   + ''.join(random.choice(letters) for i in range(4))
            n1u  = "n1u"  + ''.join(random.choice(letters) for i in range(4))
            n2   = "n2"   + ''.join(random.choice(letters) for i in range(4))
            n2u  = "n2u"  + ''.join(random.choice(letters) for i in range(4))
            m1   = "m1"   + ''.join(random.choice(letters) for i in range(4))
            m1u  = "m1u"  + ''.join(random.choice(letters) for i in range(4))
            cx   = "cx"   + ''.join(random.choice(letters) for i in range(4))
            cy   = "cy"   + ''.join(random.choice(letters) for i in range(4))
            expr  = "{phi} = atan2({cy} , {cx}) ; ".format(phi = phi, cx = cx, cy = cy)
            expr += self._dotExpression(cx, "x"+n1u,"y"+n1u,"z"+n1u,"x"+n2u,"y"+n2u,"z"+n2u) + " ; "
            expr += self._dotExpression(cy, "x"+m1u,"y"+m1u,"z"+m1u,"x"+n2u,"y"+n2u,"z"+n2u) + " ; "
            expr += self._unitvExpression("x"+m1u,"y"+m1u,"z"+m1u, "x"+m1,  "y"+m1,  "z"+m1) + " ; "
            expr += self._crossExpression("x"+m1, "y"+m1, "z"+m1,  "x"+v2u, "y"+v2u, "z"+v2u, "x"+n1u, "y"+n1u, "z"+n1u) + " ; "
            expr += self._unitvExpression("x"+n1u,"y"+n1u,"z"+n1u, "x"+n1,  "y"+n1,  "z"+n1) + " ; "
            expr += self._unitvExpression("x"+n2u,"y"+n2u,"z"+n2u, "x"+n2,  "y"+n2,  "z"+n2) + " ; "
            expr += self._unitvExpression("x"+v2u,"y"+v2u,"z"+v2u, "x"+v2,  "y"+v2,  "z"+v2) + " ; "
            expr += self._crossExpression("x"+n2, "y"+n2, "z"+n2,  "x"+v2,  "y"+v2,  "z"+v2, "x"+v3, "y"+v3, "z"+v3) + " ; "
            expr += self._crossExpression("x"+n1, "y"+n1, "z"+n1,  "x"+v1,  "y"+v1,  "z"+v1, "x"+v2, "y"+v2, "z"+v2) + " ; "
            expr += self._diffvExpression("x"+v1, "y"+v1, "z"+v1, x2, y2, z2, x1, y1, z1) + " ; "
            expr += self._diffvExpression("x"+v2, "y"+v2, "z"+v2, x3, y3, z3, x2, y2, z2) + " ; "
            expr += self._diffvExpression("x"+v3, "y"+v3, "z"+v3, x4, y4, z4, x3, y3, z3)
            return expr

    # (kf/2) (phi - phi0)^2 with a tolerance and periodicity
    def addTorsionalRestraintForce(self, particles, kphi, phi0, phitol):
        """Add a torsional restraint to the system.

        The restraint is a flat-bottom quadratic potential defined by
        a center, force constant, and tolerance. The potential has 2pi periodicity.

        Parameters
        ----------
        particles: list of ints
            The indexes of the four particles that define the torsional angle
        kphi : Quantity in (energy/angle^2) units
            Force constant of the quadratic potential
        phi0 : Quantity in angle units
            The center of the flat-bottom potential
        phitol : Quantity in angle units
            Tolerance of the flat-bottom potential

        Returns
        -------
        The reference to the CustomTorsionForce added to the system

        """
        phiforce = None
        numtorsions = 0

        expr  = "0.5*kf*( step(dm0)*max(0,db0)^2 + step(-dm0)*max(0,-da0)^2 ) ; "
        expr += self._wrapExpression("db0", "xb0", "twopi") + " ; xb0 = theta - b0 ; "
        expr += self._wrapExpression("da0", "xa0", "twopi") + " ; xa0 = theta - a0 ; "
        expr += self._wrapExpression("dm0", "xm0", "twopi") + " ; xm0 = theta - mid0 ; mid0 = (a0 + b0)/2 ; "
        expr += "twopi = 2*pi ;"
        expr += "pi = %f" % math.pi

        if self.TorsionalRestraintForce == None:
            self.TorsionalRestraintForce = mm.CustomTorsionForce(expr)
            phiforce = self.TorsionalRestraintForce
            phiforce.addPerTorsionParameter("kf")
            phiforce.addPerTorsionParameter("a0")
            phiforce.addPerTorsionParameter("b0")
            self.system.addForce(phiforce)
            numtorsions = 0
        else:
            phiforce = self.TorsionalRestraintForce
            numtorsions = phiforce.getNumTorsions()

        kf = kphi/(kilojoule_per_mole/radians**2)
        a0 = (phi0-phitol)/radians
        b0 = (phi0+phitol)/radians
        phiforce.addTorsion(particles[0], particles[1], particles[2], particles[3], np.array([kf, a0, b0], dtype=np.double))

        return phiforce

    # ** UNTESTED **
    # Adds Boresch-style Vsite restraints [J. Phys. Chem. B, Vol. 107, No. 35, 2003]
    # between 3 reference atoms of the ligand and
    # 3 reference atoms of the receptor with flat-bottom potentials
    # Reference particles of the receptor: a, b, c
    # Reference particles of the ligand:   A, B, C
    # rA = distance(A, a)
    # thetaA = angle(b, a, A)
    # phiA = dihedral(c, b, a, A)
    # (rA, thetaA, phiA are the spherical polar coordinates of A wrt the reference frame of the recepotor)
    # thetaB = angle(a, A, B)
    # phiB = dihedral(b, a, A, B)
    # phiC = dihedral(a, A, B, C)
    def _addVsiteRestraintForceBoresch(self,
                                      lig_ref_particles, rcpt_ref_particles,
                                      kfrA,  rA0,  rAtol,
                                      kfthA, thA0, thAtol,
                                      kfphA, phA0, phAtol,
                                      kfthB, thB0, thBtol,
                                      kfphB, phB0, phBtol,
                                      kfphC, phC0, phCtol):
        assert len( lig_ref_particles) == 3
        assert len(rcpt_ref_particles) == 3

        bondforce = None
        if kfrA is not None:
            bondforce = mm.CustomBondForce("0.5*kf*( step(d)*max(0,d-tol)^2 + step(-d)*max(0,-d-tol)^2 ) ; d = r - r0")
            bondforce.addPerBondParameter("kf")
            bondforce.addPerBondParameter("r0")
            bondforce.addPerBondParameter("tol")
            p0 = rcpt_ref_particles[0]
            p1 =  lig_ref_particles[0]
            kf = kfrA/(kilojoule_per_mole/nanometer**2)
            r0 = rA0/nanometer
            tol = rAtol/nanometer
            bondforce.addBond([p0,p1], np.array([kf, r0, tol ], dtype=np.double))
            self.system.addForce(bondforce)

        angleforce = None
        if kfthA is not None or kfthB is not None:
            expr  = "0.5*kf*( step(dm0)*max(0,db0)^2 + step(-dm0)*max(0,-da0)^2 ) ; "
            expr += "db0 = xb0 - pi*floor(xb0/pi + 0.5)  ; xb0 = theta - b0 ; "
            expr += "da0 = xa0 - pi*floor(xa0/pi + 0.5)  ; xa0 = theta - a0 ; "
            expr += "dm0 = xm0 - pi*floor(xm0/pi + 0.5)  ; xm0 = theta - mid0 ; mid0 = (a0 + b0)/2 ; "
            expr += "twopi = 2*pi ;"
            expr += "pi = %f" % math.pi
            angleforce = mm.CustomAngleForce(expr)
            angleforce.addPerAngleParameter("kf")
            angleforce.addPerAngleParameter("a0")
            angleforce.addPerAngleParameter("b0")
            if kfthA is not None:
                kf = kfthA/(kilojoule_per_mole/radians**2)
                a0 = (thA0 - thAtol)/radians
                b0 = (thA0 + thAtol)/radians
                p0 = rcpt_ref_particles[1]
                p1 = rcpt_ref_particles[0]
                p2 =  lig_ref_particles[0]
                angleforce.addAngle([p0,p1,p2], np.array([kf,a0,b0], dtype=np.double) )
            if kfthB is not None:
                kf = kfthB/(kilojoule_per_mole/radians**2)
                a0 = (thB0 - thBtol)/radians
                b0 = (thB0 + thBtol)/radians
                p0 = rcpt_ref_particles[0]
                p1 =  lig_ref_particles[0]
                p2 =  lig_ref_particles[1]
                angleforce.addAngle([p0,p1,p2], np.array([kf,a0,b0], dtype=np.double) )
            self.system.addForce(angleforce)

        torsforce = None
        if kfphA is not None or kfphB is not None or kfphC is not None:
            expr  = "0.5*kf*( step(dm0)*max(0,db0)^2 + step(-dm0)*max(0,-da0)^2 ) ; "
            expr += "db0 = xb0 - twopi*floor(xb0/twopi + 0.5)  ; xb0 = theta - b0 ; "
            expr += "da0 = xa0 - twopi*floor(xa0/twopi + 0.5)  ; xa0 = theta - a0 ; "
            expr += "dm0 = xm0 - twopi*floor(xm0/twopi + 0.5)  ; xm0 = theta - mid0 ; mid0 = (a0 + b0)/2 ; "
            expr += "twopi = 2*pi ;"
            expr += "pi = %f" % math.pi
            torsforce = mm.CustomTorsionForce(expr)
            torsforce.addPerTorsionParameter("kf")
            torsforce.addPerTorsionParameter("a0")
            torsforce.addPerTorsionParameter("b0")
            if kfphA is not None:
                kf = kfphA/(kilojoule_per_mole/radians**2)
                a0 = (phA0 - phAtol)/radians
                b0 = (phA0 + phAtol)/radians
                p0 = rcpt_ref_particles[2]
                p1 = rcpt_ref_particles[1]
                p2 = rcpt_ref_particles[0]
                p3 =  lig_ref_particles[0]
                torsforce.addTorsion([p0,p1,p2,p3], np.array([kf,a0,b0], dtype=np.double) )
            if kfphB is not None:
                kf = kfphB/(kilojoule_per_mole/radians**2)
                a0 = (phB0 - phBtol)/radians
                b0 = (phB0 + phBtol)/radians
                p0 = rcpt_ref_particles[1]
                p1 = rcpt_ref_particles[0]
                p2 =  lig_ref_particles[0]
                p3 =  lig_ref_particles[1]
                torsforce.addTorsion([p0,p1,p2,p3], np.array([kf,a0,b0], dtype=np.double ))
            if kfphC is not None:
                kf = kfphC/(kilojoule_per_mole/radians**2)
                a0 = (phC0 - phCtol)/radians
                b0 = (phC0 + phCtol)/radians
                p0 = rcpt_ref_particles[0]
                p1 =  lig_ref_particles[0]
                p2 =  lig_ref_particles[1]
                p3 =  lig_ref_particles[2]
                torsforce.addTorsion([p0,p1,p2,p3], np.array([kf,a0,b0],dtype=np.double ) )
            self.system.addForce(torsforce)

        return (bondforce, angleforce, torsforce)

    def addAlignmentForce(self,
                          liga_ref_particles = None, ligb_ref_particles = None,
                          kfdispl = 0.0 * kilocalorie_per_mole/angstrom**2,
                          ktheta = 0.0 * kilocalorie_per_mole,
                          kpsi = 0.0 * kilocalorie_per_mole,
                          offset = [0., 0., 0.] * angstrom):
        """
        Adds Forces to keep two ligands aligned based on reference frames defined by two sets of 3 reference atoms
        (a1, a2, a3) and (b1, b2, b3) for ligands a and b, respectively.

        See https://pubs.acs.org/doi/10.1021/acs.jcim.1c01129

        Adds a quadratic potential on the a1-b1 distance, minus an offset

        Adds a potential of the form  k [ 1 - cos(theta) ], theta being the angle between the "z" axes of the two ligands,
        defined by the a2-a1 and b2-b1 distances, to keep the axes of the two ligands aligned.

        A potential of the form  k [ 1 - cos(psi) ], psi being the dihdral angle between the planes a3-a1-b2 and a1-b2-b3
        to keep the "roll" of the two ligands in alignment. This potential is symmetrized with respect to the
        exchange of ligand labels a and b.

        Parameters
        ----------
        liga_ref_particles : list of 3 ints
            List of the indexes of the reference atoms the first ligand
        ligb_ref_particles : list of 3 ints
            List of the indexes of the reference atoms of the second ligand
        kfdispl : Quantity in energy/distance^2 units
            Force constant of the quadratic potential
        ktheta : Quantity in energy units
            Height of the 1-cos(theta) restraint potential for the angle restraint
        kpsi : Quantity in energy units
            Height of the 1-cos(psi) restraint potential for the dihedral angle restraint
        offset : List of 3 Quantity's in distance units
            Center of the restraint. The potential is applied to the a1-b1 distance after subtracting offset.
            It is used to keep the a1-b1 distance near offset.

        Returns
        -------
        A tuple of OpenMM Forces with the displacement, the angle, and dihedral restraint.

        """
        
        if liga_ref_particles and ligb_ref_particles:
            if not (len(liga_ref_particles) == 3 and len(ligb_ref_particles) == 3):
                raise ValueError("Invalid lists of reference atoms")

        expr  = " (kfdispl/2)*distsq  ; " #displacement potential
        expr += " distsq = (x1 - offx - x2)^2 + "
        expr += "          (y1 - offy - y2)^2 + "
        expr += "          (z1 - offz - z2)^2   " #square distance between b1 and a1 after displacing back b

        displforce =  mm.CustomCompoundBondForce(2, expr);
        self.system.addForce(displforce)
        displforce.addPerBondParameter("kfdispl")
        displforce.addPerBondParameter("offx")
        displforce.addPerBondParameter("offy")
        displforce.addPerBondParameter("offz")
        offv = offset / nanometer
        offx = offv[0]
        offy = offv[1]
        offz = offv[2]
        displforce.addBond([ligb_ref_particles[0],liga_ref_particles[0]] ,
                           np.array([kfdispl/((kilojoule_per_mole/nanometer**2)),
                                     offx, offy, offz], dtype=np.double) )

        expr = "(ktheta/2) * (1 - cost) ;" #theta restraint

        expr += "cost = xdn1*xdn2 + ydn1*ydn2 + zdn1*zdn2 ; "
        expr += "xdn1 = xd1/dn1 ; ydn1 = yd1/dn1 ; zdn1 = zd1/dn1 ;"
        expr += "dn1 = sqrt(xd1^2 + yd1^2 + zd1^2 ) ;"
        expr += "xd1 = x2 - x1 ; "
        expr += "yd1 = y2 - y1 ;"
        expr += "zd1 = z2 - z1 ;"
        expr += "xdn2 = xd2/dn2 ; ydn2 = yd2/dn2 ; zdn2 = zd2/dn2 ;"
        expr += "dn2 = sqrt(xd2^2 + yd2^2 + zd2^2 ) ;"
        expr += "xd2 = x4 - x3 ; "
        expr += "yd2 = y4 - y3 ; "
        expr += "zd2 = z4 - z3   "

        thetaforce = mm.CustomCompoundBondForce(4, expr)
        self.system.addForce(thetaforce)
        thetaforce.addPerBondParameter("ktheta")
        thetaforce.addBond([ligb_ref_particles[0],ligb_ref_particles[1],liga_ref_particles[0],liga_ref_particles[1] ] ,
                           np.array([ktheta/kilojoule_per_mole], dtype=np.double) )

        expr = "(kpsi/2) * (1 - cosp) ;" #psi restraint

        expr += "cosp = xvn*xwn + yvn*ywn + zvn*zwn ; "
        expr += "xvn = xv/v ; yvn = yv/v; zvn = zv/v ;"
        expr += "v = sqrt(xv^2 + yv^2 + zv^2 ) ;"
        expr += "xv = xd0 - dot01*xdn1 ;"
        expr += "yv = yd0 - dot01*ydn1 ;"
        expr += "zv = zd0 - dot01*zdn1 ;"
        expr += "dot01 =  xd0*xdn1 +  yd0*ydn1 +  zd0*zdn1 ;"
        expr += "xd0 = x3 - x1 ;"
        expr += "yd0 = y3 - y1 ;"
        expr += "zd0 = z3 - z1 ;"
        expr += "xwn = xw/w ; ywn = yw/w; zwn = zw/w ;"
        expr += "w = sqrt(xw^2 + yw^2 + zw^2 ) ;"
        expr += "xw = xd3 - dot31*xdn1 ;"
        expr += "yw = yd3 - dot31*ydn1 ;"
        expr += "zw = zd3 - dot31*zdn1 ;"
        expr += "dot31 =  xd3*xdn1 +  yd3*ydn1 +  zd3*zdn1 ;"
        expr += "xd3 = x5 - x4 ;"
        expr += "yd3 = y5 - y4 ;"
        expr += "zd3 = z5 - z4 ; "
        expr += "xdn1 = xd1/dn1 ; ydn1 = yd1/dn1 ; zdn1 = zd1/dn1 ;"
        expr += "dn1 = sqrt(xd1^2 + yd1^2 + zd1^2 ) ;"
        expr += "xd1 = x2 - x1 ; "
        expr += "yd1 = y2 - y1 ;"
        expr += "zd1 = z2 - z1 "

        psiforce = mm.CustomCompoundBondForce(5, expr)
        self.system.addForce(psiforce)
        psiforce.addPerBondParameter("kpsi")
        psiforce.addBond([ligb_ref_particles[0],ligb_ref_particles[1],ligb_ref_particles[2],
                          liga_ref_particles[0],liga_ref_particles[2] ] ,
                           np.array([0.5*kpsi/kilojoule_per_mole], dtype=np.double) )
        psiforce.addBond([liga_ref_particles[0],liga_ref_particles[1],liga_ref_particles[2],
                          ligb_ref_particles[0],ligb_ref_particles[2] ] ,
                         np.array([0.5*kpsi/kilojoule_per_mole], dtype=np.double) )
        return (displforce, thetaforce, psiforce )

    def addVsiteRestraintForceCMAngles(self,
                    lig_cm_groups = None, rcpt_cm_groups = None,
                    ktheta = None, theta0 = None, thetatol = None,
                    kphi = None, phi0 = None, phitol = None, 
                    kpsi = None, psi0 = None, psitol = None):

        """Add orientation restraining potentials between ligand and receptor.

        It is used in conjunction with addVsiteRestraintForceCMCM() to
        define the spatial and angular extent of the receptor binding
        site.

        The orientation of the ligand with respect to the receptor is
        defined with respect to a Cartesian reference frame fixed on
        the receptor defined by 3, non co-linear, reference points
        (r1, r2, r3) determined by the centroids of 3 given groups of
        atoms of the receptor.  r1 defines the origin of the frame,
        the distance r2-r1 defines the z-axis of the frame, and the
        r3-r2-r1 plane defines the xz plane of the frame.  Similarly,
        the frame of the ligand is defined by 3 reference points (l1,
        l2, l3) determined by the centroids of 3 groups of atoms of
        the ligand molecule.

        The orientation angles defined below are those between the
        elements of the frame of the ligand with respect to the frame
        of the receptor when the two frames are shifted to a common
        origin such that r1 and l1 are both at (0,0,0).

        The theta angle is defined by the angle between the z-axes of
        the two frames, that is the angle between r2-r1 and l2-l1
        directions.

        The phi angle is the dihedral between the r3-r2-r1 plane and
        the r2-l1-l2 plane. That is, essentially, the azimuthal
        spherical polar angle of the l2-l1 ligand axis with respect to
        frame fixed on the receptor.

        The psi angle is the dihedral between the r2-r1-l2 plane and
        the l1-l2-l3 plane. That is, the torsional angle of the r2,
        l1=r1, l2, l3 points, or the "twist" around the z-axis of the
        ligand.

        The restraint potential added for each angle above (theta,
        phi, psi) is a flat-bottom quadratic potential defined by
        a center, force constant, and tolerance. The potential for
        theta is in terms of cos(theta). The potentials for phi
        and psi have a periodicity of 2pi.

        Parameters
        ----------
        lig_cm_groups: list of list of ints, such as [ [0,1], [2], [3] ]
            Each element of this list is the set of atom indexes of the ligand whose centroids 
            define the reference points l1, l2, and l3, respectively.
        rcpt_cm_groups : list of list of ints, such as [ [4,5], [6], [7] ]
            Each element of this list is the set of atom indexes of the receptor whose centroids 
            define the reference points r1, r2, and r3, respectively.
        ktheta : Quantity in energy units
            Force constant of the quadratic potential for cos(theta)
        theta0 : Quantity in angle units
            The cosine of this angle is the center of the flat-bottom potential for cos(theta)
        thetatol : Quantity in angle units
            Defines the tolerance of the flat-bottom potential for cos(theta) in terms
            of cos(theta0+-thetatol)
        kphi : Quantity in (energy/angle^2) units
            Force constant of the quadratic potential for the phi angle
        phi0 : Quantity in angle units
            The center of the flat-bottom potential for the phi angle
        phitol : Quantity in angle units
            Tolerance of the flat-bottom potential for the phi angle
        kpsi : Quantity in (energy/angle^2) units
            Force constant of the quadratic potential for the psi angle
        psi0 : Quantity in angle units
            The center of the flat-bottom potential for the psi angle
        psitol : Quantity in angle units
            Tolerance of the flat-bottom potential for the psi angle

        Returns
        -------
        A tuple of OpenMM Forces for the theta, phi, and psi restraints, respectively.

        """

        assert len(lig_cm_groups) == 3
        assert len(rcpt_cm_groups) == 3

        #theta restraint
        # (kf/2) (cos(theta) - cos(theta0))^2 with a tolerance
        thetaforce = None
        expr = "0.5*kf*( step(dcost)*max(0,dcost-ctol)^2 + step(-dcost)*max(0,-dcost-ctol)^2 ) ; dcost = cost - cos0 ; "        
        expr += "cost = xdn1*xdn2 + ydn1*ydn2 + zdn1*zdn2 ; "
        expr += "xdn1 = xd1/dn1 ; ydn1 = yd1/dn1 ; zdn1 = zd1/dn1 ;"
        expr += "dn1 = sqrt(xd1^2 + yd1^2 + zd1^2 ) ;"
        expr += "xd1 = x2 - x1 ; "
        expr += "yd1 = y2 - y1 ;"
        expr += "zd1 = z2 - z1 ;"
        expr += "xdn2 = xd2/dn2 ; ydn2 = yd2/dn2 ; zdn2 = zd2/dn2 ;"
        expr += "dn2 = sqrt(xd2^2 + yd2^2 + zd2^2 ) ;"
        expr += "xd2 = x4 - x3 ; "
        expr += "yd2 = y4 - y3 ; "
        expr += "zd2 = z4 - z3   "
        if ktheta is not None:
            if self.CMAngleThetaForce is None:
               self.CMAngleThetaForce = mm.CustomCentroidBondForce(4, expr)
               thetaforce = self.CMAngleThetaForce
               numthetagroups = 0
               thetaforce.addPerBondParameter("kf")
               thetaforce.addPerBondParameter("cos0")
               thetaforce.addPerBondParameter("ctol")
               self.system.addForce(thetaforce)
            else:
                thetaforce = self.CMAngleThetaForce
                numthetagroups = thetaforce.getNumGroups()
            thetaforce.addGroup(rcpt_cm_groups[0])
            thetaforce.addGroup(rcpt_cm_groups[1])
            thetaforce.addGroup( lig_cm_groups[0])
            thetaforce.addGroup( lig_cm_groups[1])
            kf = ktheta/(kilojoule_per_mole)
            cos0 = math.cos(theta0/radians)
            a0 = max(0, (theta0-thetatol)/radians)
            b0 = min(math.pi, (theta0+thetatol)/radians)
            ctol = abs(math.cos(a0) - math.cos(b0))
            groupids = [ numthetagroups + i for i in range(4) ]
            thetaforce.addBond(groupids, np.array([kf, cos0, ctol], dtype=np.double) )

        #phi restraint
        # (kf/2) (phi - phi0)^2 with a tolerance and periodicity
        phiforce = None
        expr  = "0.5*kf*( step(dm0)*max(0,db0)^2 + step(-dm0)*max(0,-da0)^2 ) ; "
        expr += self._wrapExpression("db0", "xb0", "twopi") + " ; xb0 = psi - b0 ; "
        expr += self._wrapExpression("da0", "xa0", "twopi") + " ; xa0 = psi - a0 ; "
        expr += self._wrapExpression("dm0", "xm0", "twopi") + " ; xm0 = psi - mid0 ; mid0 = (a0 + b0)/2 ; "
        expr += self._dihedralExpression("phi", "xr2", "yr2", "zr2", "xr1", "yr1", "zr1", "0", "0", "0", "xl1", "yl1", "zl1") + " ; "
        expr += self._diffvExpression("xr2", "yr2", "zr2", "x3", "y3", "z3",  "x1", "y1", "z1") + " ; "
        expr += self._diffvExpression("xr1", "yr1", "zr1", "x2", "y2", "z2",  "x1", "y1", "z1") + " ; "
        expr += self._diffvExpression("xl1", "yl1", "zl1", "x5", "y5", "z5",  "x4", "y4", "z4") + " ; "
        expr += "twopi = 2*pi ;"
        expr += "pi = %f" % math.pi
        if kphi is not None:
            if self.CMAnglePhiForce is None:
               self.CMAnglePhiForce = mm.CustomCentroidBondForce(5, expr)
               phiforce = self.CMAnglePhiForce
               numphigroups = 0
               phiforce.addPerBondParameter("kf")
               phiforce.addPerBondParameter("a0")
               phiforce.addPerBondParameter("b0")
               self.system.addForce(phiforce)
            else:
                phiforce = self.CMAnglePhiForce
                numphigroups = phiforce.getNumGroups()
            phiforce.addGroup(rcpt_cm_groups[0])#1
            phiforce.addGroup(rcpt_cm_groups[1])#2
            phiforce.addGroup(rcpt_cm_groups[2])#3
            phiforce.addGroup( lig_cm_groups[0])#4
            phiforce.addGroup( lig_cm_groups[1])#5
            kf = kphi/(kilojoule_per_mole/radians**2)
            a0 = (phi0-phitol)/radians
            b0 = (phi0+phitol)/radians
            groupids = [ numphigroups + i for i in range(5) ]
            phiforce.addBond(groupids, np.array([kf, a0, b0], dtype=np.double)  )

        #psi restraint
        # (kf/2) (psi - psi0)^2 with a tolerance and periodicity
        psiforce = None
        expr  = "0.5*kf*( step(dm0)*max(0,db0)^2 + step(-dm0)*max(0,-da0)^2 ) ; "
        expr += self._wrapExpression("db0", "xb0", "twopi") + " ; xb0 = psi - b0 ; "
        expr += self._wrapExpression("da0", "xa0", "twopi") + " ; xa0 = psi - a0 ; "
        expr += self._wrapExpression("dm0", "xm0", "twopi") + " ; xm0 = psi - mid0 ; mid0 = (a0 + b0)/2 ; "
        expr += self._dihedralExpression("psi", "xr1", "yr1", "zr1", "0", "0", "0", "xl1", "yl1", "zl1", "xl2", "yl2", "zl2") + " ; "
        expr += self._diffvExpression("xr1", "yr1", "zr1", "x2", "y2", "z2",  "x1", "y1", "z1") + " ; "
        expr += self._diffvExpression("xl1", "yl1", "zl1", "x4", "y4", "z4",  "x3", "y3", "z3") + " ; "
        expr += self._diffvExpression("xl2", "yl2", "zl2", "x5", "y5", "z5",  "x3", "y3", "z3") + " ; "
        expr += "twopi = 2*pi ;"
        expr += "pi = %f" % math.pi
        if kpsi is not None:
            if self.CMAnglePsiForce is None:
               self.CMAnglePsiForce = mm.CustomCentroidBondForce(5, expr)
               psiforce = self.CMAnglePsiForce
               numpsigroups = 0
               psiforce.addPerBondParameter("kf")
               psiforce.addPerBondParameter("a0")
               psiforce.addPerBondParameter("b0")
               self.system.addForce(psiforce)
            else:
                psiforce = self.CMAnglePsiForce
                numpsigroups = psiforce.getNumGroups()
            psiforce.addGroup(rcpt_cm_groups[0])#1
            psiforce.addGroup(rcpt_cm_groups[1])#2
            psiforce.addGroup( lig_cm_groups[0])#3
            psiforce.addGroup( lig_cm_groups[1])#4
            psiforce.addGroup( lig_cm_groups[2])#5
            kf = kpsi/(kilojoule_per_mole/radians**2)
            a0 = (psi0-psitol)/radians
            b0 = (psi0+psitol)/radians
            groupids = [ numpsigroups + i for i in range(5) ]
            psiforce.addBond(groupids, np.array([kf, a0, b0], dtype=np.double) )
        return (thetaforce, phiforce, psiforce)

    def addPosRestraints(self, particles, refpos, fc = 25.0 * kilocalorie_per_mole/angstrom**2, tol = 0.5 * angstrom, periodic = True):
        """
        Applies flat-bottom, quadratic position restraints to a set of
        atoms based on a set of reference positions.

        Parameters
        ----------
        particles: list of ints
            List of the indexes of the atoms to restrain
        refpos : list of Vec3
            Centers of the flat-bottom restraints. These are reference 
            positions for all atoms of the system, 
            not just those that are restrained.
        fc : Quantitiy in (energy/distance^2) units
            Force constant of the quadratic flat-bottom potential
        tol : Quantity in distance units
            Tolerance of the quadratic flat-bottom potential

        Returns
        -------
        The Force added to the System

        """

        if not particles or len(particles) == 0:
            return None
        if periodic:
            posrestforce = mm.CustomExternalForce("0.5*fc*select(step(dist-tol), (dist-tol)^2, 0); dist = periodicdistance(x,y,z,x0,y0,z0)")
        else:
            posrestforce = mm.CustomExternalForce("0.5*fc*select(step(dist-tol), (dist-tol)^2, 0); dist = sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)")

        posrestforce.addPerParticleParameter("x0")
        posrestforce.addPerParticleParameter("y0")
        posrestforce.addPerParticleParameter("z0")
        posrestforce.addPerParticleParameter("fc")
        posrestforce.addPerParticleParameter("tol")

        self.system.addForce(posrestforce)

        for p in particles:
            x0 = refpos[p][0]/nanometer
            y0 = refpos[p][1]/nanometer
            z0 = refpos[p][2]/nanometer
            fc1 = fc/(kilojoule_per_mole/nanometer**2)
            tol1 = tol/nanometer
            posrestforce.addParticle(p, np.array([x0, y0, z0, fc1, tol1], dtype=np.double)  )
        return posrestforce

    #fix zero LJs
    def fixZeroLJParams(self, force, minsigma = 0.1 * angstrom, minepsilon = 1.e-4 * kilocalories_per_mole):
        """
        Set the Lennard-Jones parameters to some minimum values

        This is used for unphysical force fields that do not have
        a ground state because of "naked" charge centers that can
        sit on other charge centers without penalty.

        Parameters
        ----------
        force : OpenMM force
            Non-bonded Force to operate upon, other Forces are skipped
        minsigma : Quantity in distance units
            Minimum value of the sigma LJ parameter to assign
        minepsilone : Quantity in energy units
            Minimum value of the epsilon LJ parameter to assign

        """

        small = 1.0e-6
        if isinstance(force, mm.NonbondedForce):
            for i in range(force.getNumParticles()):
                nbparam = force.getParticleParameters(i)
                charge = nbparam[0]
                sigmaLJ = nbparam[1]
                epsilonLJ = nbparam[2]
                if sigmaLJ._value < small and epsilonLJ._value < small:
                    sigmaLJ = minsigma
                    epsilonLJ = minepsilon
                    force.setParticleParameters(i, charge, sigmaLJ, epsilonLJ)

    def softCorePertE(self, u,  umax, ub, a):
        usc = u
        if u > ub :
            gu = (u-ub)/(a*(umax-ub)) #this is y/alpha
            zeta = 1. + 2.*gu*(gu + 1.)
            zetap = np.power( zeta , a)
            usc = (umax-ub)*(zetap - 1.)/(zetap + 1.) + ub
        return usc

def separate14(nbforce):
    ONE_4PI_EPS0 = 138.93545764438198
    expr = "4*epsilon*((sigma/r)^12 - (sigma/r)^6) + ONE_4PI_EPS0*chargeprod/r;" + \
        "ONE_4PI_EPS0 = {:f};".format(ONE_4PI_EPS0)
    force14 = mm.CustomBondForce(expr)
    force14.addPerBondParameter('chargeprod')
    force14.addPerBondParameter('sigma')
    force14.addPerBondParameter('epsilon')
    for i in range(nbforce.getNumExceptions()):
        p1, p2, chargeprod, sigma, epsilon = nbforce.getExceptionParameters(i)
        if(math.fabs(chargeprod._value) > 1.e-8 or math.fabs(epsilon._value) > 1.e-8):
            force14.addBond(p1, p2, [chargeprod, sigma, epsilon])
            nbforce.setExceptionParameters(i, p1, p2, 0.0, 0.1, 0.0)#turn off 1-4 term
    return force14

def residue_is_solvent(res): # called in abfe/rbfe_structprep.py
    ion_resnames = ['POT','SOD','CLA','NA+','K+','CL-','F-','CA','MG','CL','NA','K','F']
    wat_resnames = ['HOH','TIP3','WAT','TIP4','OPC','TIP5'] # sometimes we have 'TIP3V'      
    if res.name.upper() in (ion_resnames + wat_resnames):
        return True
    return False

