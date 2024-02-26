# setup.py
# Install script of ASyncRE module
# Copyright (C) 2021 Emilio Gallicchio
# E-mail: emilio.gallicchio@gmail.com
#
# This software is licensed under the terms of the GNU General Public License
# http://opensource.org/licenses/GPL-3.0
#
#    This program is free software: you can redistribute it and/or
#    modify it under the terms of the GNU General Public License
#    version 3 as published by the Free Software Foundation.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

from setuptools import setup
from async_re import __version__ as VERSION

NAME = 'atom-openmm'

MODULES = 'async_re', 'ommreplica', 'ommsystem', 'ommworker', 'local_openmm_transport', 'transport', 'gibbs_sampling', 'openmm_async_re'


SCRIPTS = 'abfe_explicit.py', 'rbfe_explicit.py', 'abfe_structprep.py', 'rbfe_structprep.py','abfe_explicit_sync.py', 'rbfe_explicit_sync.py', 'make_atm_system_from_pdb.py', 'make_atm_system_from_Amber.py', 'make_atm_system_from_rcpt_lig.py'

REQUIRES = 'configobj', 'numpy' 

DESCRIPTION = 'Asynchronous Replica Exchange with OpenMM'

AUTHOR = 'Emilio Gallicchio'

AUTHOR_EMAIL = 'emilio.gallicchio@gmail.com'

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      py_modules=MODULES,
      scripts=SCRIPTS,
      packages=['sync', 'utils'],
      package_data={'utils': ['logging.conf']},
      requires=REQUIRES
     )
