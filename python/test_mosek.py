# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 16:17:06 2018

@author: gualandi
"""

from mosek import *;

x = [ 0.0 ]

with Env() as env:                            # Create Environment
  with env.Task(0, 1) as task:                # Create Task
    task.appendvars(1)                          # 1 variable x
    task.putcj(0, 1.0)                          # c_0 = 1.0
    task.putvarbound(0, boundkey.ra, 2.0, 3.0)  # 2.0 <= x <= 3.0
    task.putobjsense(objsense.minimize)         # minimize

    task.optimize()                      # Optimize

    task.getxx(soltype.itr, x)                  # Get solution
    print("Solution x = {}".format(x[0]))       # Print solution