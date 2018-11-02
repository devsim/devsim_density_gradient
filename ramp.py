# Copyright 2013 Devsim LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys
from devsim import *
import devsim
#from simple_physics import *

def rampparam(device, region, param, stop, step_size, min_step, solver_params=None, callback=None):
  '''
    Ramps param with assignable callback function
  '''
  start_param=get_parameter(device=device, region=region, name=param)
  if (start_param < stop):
    step_sign=1
  else:
    step_sign=-1
  last_param=start_param
  while(abs(last_param - stop) > min_step):
    print("last end %e %e" % (last_param, stop))
    next_param=last_param + step_sign * step_size
    if next_param < stop:
      next_step_sign=1
    else:
      next_step_sign=-1

    if next_step_sign != step_sign:
      next_param=stop
      print("setting to last param %e" % (stop))
      print("setting next param %e" % (next_param))
    set_parameter(device=device, region=region, name=param, value=next_param)
    try:
      solve(**solver_params)
    except devsim.error as msg:
      if str(msg).find("Convergence failure") != 0:
        raise
      set_parameter(device=device, region=region, name=param, value=last_param)
      step_size *= 0.5
      print("setting new step size %e" % (step_size))
      #raise NameError("STOP")
      if step_size < min_step:
        raise RuntimeError("Min step size too small")
      continue
    last_param=next_param
    if callback:
      callback()

def rampbias(device, contact, end_bias, step_size, min_step, max_iter, rel_error, abs_error, callback):
  '''
    Ramps bias with assignable callback function
  '''
  start_bias=get_parameter(device=device, name=GetContactBiasName(contact))
  if (start_bias < end_bias):
    step_sign=1
  else:
    step_sign=-1
  last_bias=start_bias
  while(abs(last_bias - end_bias) > min_step):
    print(("last end %e %e") % (last_bias, end_bias))
    next_bias=last_bias + step_sign * step_size
    if next_bias < end_bias:
      next_step_sign=1
    else:
      next_step_sign=-1

    if next_step_sign != step_sign:
      next_bias=end_bias
      print("setting to last bias %e" % (end_bias))
      print("setting next bias %e" % (next_bias))
    set_parameter(device=device, name=GetContactBiasName(contact), value=next_bias)
    try:
      solve(type="dc", absolute_error=abs_error, relative_error=rel_error, maximum_iterations=max_iter)
    except devsim.error as msg:
      if msg[0].find("Convergence failure") != 0:
        raise
      set_parameter(device=device, name=GetContactBiasName(contact), value=last_bias)
      step_size *= 0.5
      print("setting new step size %e" % (step_size))
      if step_size < min_step:
        raise RuntimeError("Min step size too small")
      continue
    print("Succeeded")
    last_bias=next_bias
    callback(device)

def printAllCurrents(device):
  '''
    Prints all contact currents on device
  '''
  for c in get_contact_list(device=device):
    PrintCurrents(device, c)

