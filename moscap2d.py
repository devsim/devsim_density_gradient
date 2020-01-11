from devsim import *
def instantiate(filename, device, interface_siox, region_ox, region_si, contact_ox, contact_si, ox_material, si_material):
    create_gmsh_mesh (mesh=device, file=filename)
    add_gmsh_region  (mesh=device, gmsh_name="oxide",   region=region_ox, material=ox_material)
    add_gmsh_region  (mesh=device, gmsh_name="silicon", region=region_si, material=si_material)
    add_gmsh_contact (mesh=device, gmsh_name="gate",    region=region_ox, material="metal", name=contact_ox)
    add_gmsh_contact (mesh=device, gmsh_name="sub",     region=region_si, material="metal", name=contact_si)
    add_gmsh_interface (mesh=device, gmsh_name="interface", region0=region_ox, region1=region_si, name=interface_siox)
    finalize_mesh    (mesh=device)
    create_device    (mesh=device, device=device)

def run():
    test_opts = {
        "filename" : "moscap2d.msh",
      "device" : "MyDevice",
      "region_ox" : "MyOxRegion",
      "region_si" : "MySiRegion",
      "interface_siox" : "MySiOx",
      "contact_ox" : "top",
      "contact_si" : "bot",
      "ox_material" : "Ox",
      "si_material" : "Si"
    }
    instantiate(**test_opts)
    write_devices(file="moscap2d_devsim.msh", type="devsim")
    write_devices(file="moscap2d.tec", type="tecplot")

#write_devices(file="moscap.dsm", type="devsim")
if __name__ == '__main__':
    run()

