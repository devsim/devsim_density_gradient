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

####
#### cap2.py
#### tests physics of cap made of two insulating regions
####
from devsim import *
def make_moscap(device, interface_siox, region_ox, region_si, contact_ox, contact_si, tox, tsi, ox_material, si_material):
  '''
    need to generalize, but good enough for our examples
  '''
  create_1d_mesh(mesh="moscap")
  if False:
    add_1d_mesh_line(mesh="moscap", pos=tox,   ps=1e-7, tag="d0")
    add_1d_mesh_line(mesh="moscap", pos=0.0, ns=1e-9, ps=1e-7, tag="d1")
    if tsi > 1e-6:
      add_1d_mesh_line(mesh="moscap", pos=-1e-6,   ps=5e-9)
    add_1d_mesh_line(mesh="moscap", pos=-tsi,   ps=1e-7, tag="d2")
    add_1d_contact  (mesh="moscap", name=contact_ox,    tag="d0", material="metal")
    add_1d_contact  (mesh="moscap", name=contact_si,    tag="d2", material="metal")
    add_1d_interface(mesh="moscap", name=interface_siox, tag="d1")
    add_1d_region   (mesh="moscap", material=ox_material, region=region_ox, tag1="d0", tag2="d1")
    add_1d_region   (mesh="moscap", material=si_material, region=region_si, tag1="d1", tag2="d2")
  else:
    add_1d_mesh_line(mesh="moscap", pos=-tox,   ps=1e-7, tag="d0")
    add_1d_mesh_line(mesh="moscap", pos=0.0, ps=1e-9, ns=1e-7, tag="d1")
    if tsi > 1e-6:
      add_1d_mesh_line(mesh="moscap", pos=1e-6,   ps=5e-9)
    add_1d_mesh_line(mesh="moscap", pos=tsi,   ps=1e-7, tag="d2")
    add_1d_contact  (mesh="moscap", name=contact_ox,    tag="d0", material="metal")
    add_1d_contact  (mesh="moscap", name=contact_si,    tag="d2", material="metal")
    add_1d_interface(mesh="moscap", name=interface_siox, tag="d1")
    add_1d_region   (mesh="moscap", material=ox_material, region=region_ox, tag1="d0", tag2="d1")
    add_1d_region   (mesh="moscap", material=si_material, region=region_si, tag1="d1", tag2="d2")
  finalize_mesh(mesh="moscap")
  create_device(mesh="moscap", device=device)

#create_contact_from_interface(name="mid", material="metal", device=device, region="MySiRegion", interface=interface)
#print [x for x in get_node_model_values(device=device, region="MySiRegion", name="AtContactNode") if x > 0]
#print [x for x in get_node_model_values(device=device, region="MyOxRegion", name="AtContactNode") if x > 0]
#print [x for x in get_node_model_values(device=device, region="MyOxRegion", name="SurfaceArea") if x > 0]
#print [x for x in get_node_model_values(device=device, region="MySiRegion", name="SurfaceArea") if x > 0]

#write_devices(file="moscap.dsm", type="devsim")

