#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 10:41:17 2020

@author: Elena Zamaraeva
"""

import numpy as np
import json

"""
    This file contains a state class used in RL
"""


class State:

    def __init__(self, energy, struct=[], old_energy=0, unique=True):
        """
            A state has the following features:
            energy: the energy of the current structure
            old_energy: the energy of the previous structure (can be used for zero-reward case)
            isunique: indicates if the energy was not met before
            struct: the structure desccription, usually the array of the atoms and their positions
            depending on the representation in the main CSP code
        """
        self.energy = energy
        self.struct = struct
        self.old_energy = old_energy
        self.isunique = unique

    def __str__(self):
        return str(self.__class__) + ": " + str(self.__dict__)

    def sl(self):
        """
            Returns the stack length for MC-EMMA structures
        """
        return len(self.struct)

    def struct_to_str(self):

        if type(self.struct) is np.ndarray:
            struct = self.struct.tolist()
        else:
            struct = self.struct

        return json.dumps(struct)

    def copy(self):

        return State(self.energy, self.struct.copy(), old_energy=self.old_energy, unique=self.isunique)
