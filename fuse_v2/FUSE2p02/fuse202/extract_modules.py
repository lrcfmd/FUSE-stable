#!/usr/bin/env python
# coding: utf-8

# In[226]:

from ase.io import *
from ase.visualize import *
from ase import Atoms
import os
import sys
import glob
from fuse202 import bond_table
import numpy

bondtable=numpy.load("bondtable.npz",allow_pickle=True)
bondtable=bondtable['bond_table'].item()

def extract_module(input_files,bondtable,man_ap=None,max_z=0.8,max_xy=0.8):
        try:
            atoms=[]
            
            #0. Pre-process all of the cif files to make sure that all co-ordinates are between 0 and 1
            
            for i in input_files:
            
                at=read(i)
            
                at.set_scaled_positions(at.get_scaled_positions() )
            
                atoms.append(at)
            
            full_structure=atoms[0].copy()
            
            modules=[]
            
            for i in atoms:
            	  
                # 1. Work out the unit cell dimensions for each structure in terms of sub modules, using atomic radii, or through a manual sizing                   
                
                form=i.get_chemical_symbols()
                if man_ap == None:
                    initial_guess=0
                    for j in form:
                        initial_guess+=bondtable[j][0][2] 
                    initial_guess=initial_guess/(len(i)/4)
                    initial_guess=round(initial_guess,6)
                    
                else:
                    initial_guess=round(man_ap,6)
                            
                sub_modules=[0,0,0,90,90,90]
                cell=i.cell.cellpar()
                
                sub_modules=[round(cell[0]/initial_guess),round(cell[1]/initial_guess),round((cell[2]/initial_guess)*2),cell[3],cell[4],cell[5]]
                
                if sub_modules[0]==0:
                	 sub_modules[0]=1
                
                if sub_modules[1]==0:
                	 sub_modules[1]=1
                
                if sub_modules[2]==0:
                	 sub_modules[1]=1
            
                sub_module_cell=[round(initial_guess,6),round(initial_guess,6),round(initial_guess/2,6),round(cell[3]),round(cell[4]),round(cell[5])]
                
                 
                # 2. Now that we have our unit cell, break each structure down into it's component sub-modules
                
                n_sub=sub_modules[0]*sub_modules[1]*sub_modules[2]
                n_sub_made=0        
                position=[0,0,0] # position in sub-modules, rows alon x, then y, then layers along z
                
                #will break down structure, starting at (0,0,0) moving in rows along y, then up through z.
            
                co_ords=[ [ (  cell[0]/sub_modules[0]  * (position[0]/cell[0]) ),   (cell[0] / sub_modules[0]) * ((position[0]+1)/cell[0]) ], 
                          [ (  cell[1]/sub_modules[1]  * (position[1]/cell[1]) ),   (cell[1] / sub_modules[1]) * ((position[1]+1)/cell[1]) ], 
                          [ (  cell[2]/sub_modules[2]  * (position[2]/cell[2]) ),   (cell[2] / sub_modules[2]) * ((position[2]+1)/cell[2]) ] ]
                
                for j in range(n_sub):
                    
                    # After the first sub module, assign the new set of co-ordinates to pick out for sub-modules:
            
                    co_ords=[ [ (  cell[0]/sub_modules[0]  * (position[0]/cell[0]) ) ,  (cell[0] / sub_modules[0]) * ((position[0]+1)/cell[0]) ], 
                              [ (  cell[1]/sub_modules[1]  * (position[1]/cell[1]) ) ,  (cell[1] / sub_modules[1]) * ((position[1]+1)/cell[1]) ], 
                              [ (  cell[2]/sub_modules[2]  * (position[2]/cell[2]) ) ,  (cell[2] / sub_modules[2]) * ((position[2]+1)/cell[2]) ] ]
                                        
                    #sub_module_cell
                    mod=Atoms()
                    mod.pbc=['True','True','True']
                    mod.cell=sub_module_cell
                    
                    for k in i:
                        #print(k)
                        sp = k.scaled_position
                        if sp[0] >= co_ords[0][0]:
                            if sp[0] < co_ords[0][1]:
                                if sp[1] >= co_ords[1][0]:
                                    if sp[1] < co_ords[1][1]:
                                        if sp[2] >= co_ords[2][0]:
                                            if sp[2] < co_ords[2][1]:
                                                # before we put it into the sub module, need to work out what the co-ordiates should be relative to the origin of the sub-mod
                                                zero_point=[cell[0]*co_ords[0][0], cell[1]*co_ords[1][0], cell[2]*co_ords[2][0] ]
                                                pos=k.position-zero_point
                                                # work out the new position as a fraction of the co-ordinate range we're looking at
                                                rng=[(co_ords[0][1] - co_ords[0][0])*cell[0] , (co_ords[1][1] - co_ords[1][0])*cell[1] , (co_ords[2][1] - co_ords[2][0])*cell[2] ] 
                                                pos=[pos[0]/rng[0], pos[1]/rng[1], pos[2]/rng[2]]
                                                new_pos=[pos[0]*sub_module_cell[0], pos[1]*sub_module_cell[1], pos[2]*sub_module_cell[2]]
                                                k.position=new_pos
                                                
                                                mod=mod + k
                                                
                    
                    # 3. need to reset the co-ordinates within the module to be within the sub-module
                    
                    pos=mod.get_scaled_positions(wrap=True)
                    mod.set_scaled_positions(mod.get_scaled_positions(wrap=True))
                    modules.append(mod)
                    n_sub_made+=1
                    
                    # 4. Update the position variable to move to the next sub-module to extract
                    
                    if j != n_sub-1:
                        
                        # see if we need to continue along x
                        if position[0] < sub_modules[0]-1: 
                            position[0]+=1
                            continue
                            
                        if position[0] == sub_modules[0] -1:
                            # see if we need to continue along y, or move up a row
                            
                            if position[1] == sub_modules[1]-1:
                                position[2]+= 1
                                position[0]=0
                                position[1]=0
                                continue
                                
                            # if not, continue along the layer
                            
                            if position[1] < sub_modules[1]-1:
                                position[1] += 1
                                position[0] = 0
                             
                        
                            
                        # Need to move the range of co-ordinates that we're searching for for the next sub-module
                        
                               
            structure={'modules':modules,'sub module cell':sub_module_cell,'shape in submods':sub_modules,'nmods':n_sub,'ap':initial_guess,'atoms':full_structure}
            
            #if there's a specified limit on the maximum fractional co-ordinate for a sub module, go through and scale the sub_module positions:
            scale=False
            if max_z != None:
                scale=True
            	
            if max_xy!= None:
                scale=True
            	
            if scale == True:
                mods=structure['modules'].copy()
                s=[1,1,1]
                if max_z != None:
                	s[2]=max_z
                if max_xy != None:
                	s[0]=max_xy
                	s[1]=max_xy    
                	
                for i in list(range(len(mods))):
                	temp=mods[i].copy()
                	#print(temp)
                	#print(temp.positions)
                	for j in list(range(len(temp))):
                		temp[j].position=[temp[j].position[0]*s[0],temp[j].position[1]*s[1],temp[j].position[2]*s[2]]
                	#print(temp.positions)
                	
                	mods[i]=temp.copy()
                	
                	#sys.exit()
                structure['modules']=mods.copy()	
            
            
            return structure

        except:
        	  structure=None
        	  return structure


