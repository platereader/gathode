"""
Parser classes to read exported files by different platereaders, used by GATHODE.

Growth Analysis Tool for High-throughput Optical Density Experiments
(GATHODE) Plate class.
"""

# GATHODE  Growth Analysis Tool
#          for High-throughput Optical Density Experiments
#
# Copyright (C) 2014 Nils Christian
#
# This file is part of GATHODE.
#
# GATHODE is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# GATHODE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with GATHODE.  If not, see <http://www.gnu.org/licenses/>.

__all__ = ['tecan','bioscreen']

import sys
import inspect

from . import *

def getModulesOfNamespace(namespace,handledModules=[],depth=0,orignamespace=None,exclude=None):
    if orignamespace is None:
        orignamespace = namespace.__name__
    if exclude is None:
        exclude = []
    modules=set()
    for name, obj in inspect.getmembers(namespace):
        if inspect.ismodule(obj) and str(obj) not in modules and name not in exclude and obj.__name__.startswith(str(orignamespace)):
            modules.update(set([obj]))
            newclasses = getModulesOfNamespace(obj,handledModules,depth=depth+1,orignamespace=orignamespace)
            modules.update(newclasses)
    return modules

def modulenameToModule(modules,replace=None,replacewith='',lower=False):
    name2module = {}
    for md in modules:
        name = md.__name__
        if lower:
            name = name.lower()
        if replace is not None:
            name = name.replace(replace,replacewith)
        name2module[name]=md
    return name2module

parser2module = modulenameToModule(
    list(getModulesOfNamespace(sys.modules[__name__],exclude=['utils'])),
    replace='platereader.parser.',
    lower=True)
