import devsim
from dg_physics import *
from dg_common import *
import moscap
import numpy
import ramp
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import sys

def latex_float(f):
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        float_str = ''
        if base == '-1':
          float_str = r"-$10^{{{0}}}$".format(int(exponent))
        elif base == '1':
          float_str = r"$10^{{{0}}}$".format(int(exponent))
        else:
          float_str = r"${0} \cdot 10^{{{1}}}$".format(base, int(exponent))
    return float_str

def format_doping(f):
  s = latex_float(f)
  if s[0] == '-':
    s = r'N$_A$=' + s[1:]
  else:
    s = r'N$_D$=' + s
  s = s + r' (#/cm$^3$)'
  return s


def simulate_charge(device, contact, equation, solver_params):
  factor=1e7 #from F/cm^2 to fF/um^2
  dv = 0.001
  v1 = get_parameter(device=device, name=GetContactBiasName(contact))
  q1 = get_contact_charge(device=device, contact=contact, equation="PotentialEquation")
  v2 = v1 + 0.001
  set_parameter(name=GetContactBiasName(contact), value=v2)
  devsim.solve(**solver_params)
  q2 = get_contact_charge(device=device, contact=contact, equation="PotentialEquation")
  return (v1, (factor*(q2-q1)/dv))



#assume body nodes have two edges
def debug_plot(device, region):
  nodes = numpy.array(get_node_model_values(device=device, region=region, name='node_index'))
  node_volume = numpy.array(get_node_model_values(device=device, region=region, name='NodeVolume'))
  num_nodes = len(nodes)
  flux_term     = numpy.zeros(num_nodes)
  edge_vol_term = numpy.zeros(num_nodes)
  x             = abs(numpy.array(get_node_model_values(device=device, region=region, name="x")))*1e7
  node_vol_term = numpy.array(get_node_model_values(device=device, region=region, name="volume_term"))
  pos_edge = []
  neg_edge = []
  for i in range(num_nodes):
    pos_edge.append([])
    neg_edge.append([])

  # edge based quantities
  edge_from_node_model(device=device, region=region, node_model='node_index')
  n0_indexes = get_edge_model_values(device=device, region=region, name='node_index@n0')
  n1_indexes = get_edge_model_values(device=device, region=region, name='node_index@n1')
  # need to get sign based on index
  num_edges = len(n0_indexes)
  zero_earray = [0.0]*num_edges

  # populate edges
  for ei in range(num_edges):
    pos_edge[int(n0_indexes[ei])].append(ei)
    neg_edge[int(n1_indexes[ei])].append(ei)

  ft = get_edge_model_values(device=device, region=region, name="surface_term")
  st = get_edge_model_values(device=device, region=region, name="edge_volume_term")
  for ni in range(num_nodes):
    # don't want a boundary yet
    for v in pos_edge[ni]:
      flux_term[ni]     += ft[v]
      edge_vol_term[ni] += st[v]
    for v in neg_edge[ni]:
      flux_term[ni]     -= ft[v]
      edge_vol_term[ni] += st[v]
  node_vol_term = node_vol_term / node_volume
  edge_vol_term = edge_vol_term / node_volume
  flux_term     = flux_term / node_volume
  net_sum = flux_term + edge_vol_term+ node_vol_term
  fig=plt.figure()
  plt.plot(x, node_vol_term, label="node vol")
  plt.plot(x, edge_vol_term, label="edge vol")
  plt.plot(x, flux_term, label="flux")
  plt.plot(x, net_sum, 'k+-', label="sum")
  plt.legend(loc='best')
  plt.ylim(-1, 1.0)
  plt.xlim(0, 20)
  plt.xlabel(xstring)
  plt.ylabel('eq terms (eV)')
  return fig

def setup_dd_si(device, region):
  # this is our solution variable
  CreateSolution(device, region, "Potential")
  CreateSolution(device, region, "Electrons")
  CreateSolution(device, region, "Holes")


  #These are needed for the Arora model
  CreateBandEdgeModels(device, region, ("Potential","Le","Lh"))

  #these are needed for velocity saturation
  #CreateEField(device, region)
  #CreateDField(device, region)
  opts = CreateAroraMobilityLF(device, region)
  #opts = CreateHFMobility(device, region, **opts)

  CreateSiliconDriftDiffusion(device, region, **opts)
  for i in get_contact_list(device=device):
    r=get_region_list(device=device, contact=i)[0]
    if r == region:
      set_parameter(name=GetContactBiasName(i), value=0.0)
      CreateSiliconDriftDiffusionContact(device, region, i, opts['Jn'], opts['Jp'])
  return opts

def ActivateLe():
  set_parameter(device=device, region=region_si, name="Gamman", value=3)
  set_parameter(device=device, region=region_ox, name="Gamman", value=1)
  #setup_dg_equation(device, region_si, 'Lambda_e', 'Le', '', '', 'Le_eqn')
  #setup_dg_equation(device, region_ox, 'Lambda_e', 'Le', '', '', 'Le_eqn')
  #setup_dg_equation(device, region_si, 'Lambda_e', 'Le', 'del_log_n1', '', 'Le_eqn')
  #setup_dg_equation(device, region_ox, 'Lambda_e', 'Le', 'del_log_n1', '', 'Le_eqn')
  setup_dg_equation(device, region_si, 'Lambda_e', 'Le', 'del_log_n1', 'del_log_n2', 'Le_eqn')
  setup_dg_equation(device, region_ox, 'Lambda_e', 'Le', 'del_log_n1', 'del_log_n2', 'Le_eqn')
  for c in ('top', 'bot'):
    setup_dg_contact(device, c, 'Lambda_e', 'Le') 
  setup_dg_interface(device, interface_siox, 'Lambda_e', 'Le')


def ActivateLh():
  set_parameter(device=device, region=region_si, name="Gammap", value=3)
  set_parameter(device=device, region=region_ox, name="Gammap", value=1)
  #setup_dg_equation(device, region_si, 'Lambda_h', 'Lh', 'del_log_p1', '', 'Lh_eqn')
  #setup_dg_equation(device, region_ox, 'Lambda_h', 'Lh', 'del_log_p1', '', 'Lh_eqn')
  setup_dg_equation(device, region_si, 'Lambda_h', 'Lh', 'del_log_p1', 'del_log_p2', 'Lh_eqn')
  setup_dg_equation(device, region_ox, 'Lambda_h', 'Lh', 'del_log_p1', 'del_log_p2', 'Lh_eqn')
  for c in ('top', 'bot'):
    setup_dg_contact(device, c, 'Lambda_h', 'Lh') 
  setup_dg_interface(device, interface_siox, 'Lambda_h', 'Lh')

def ActivateLe2():
  set_parameter(device=device, region=region_si, name="Gamman", value=3)
  set_parameter(device=device, region=region_ox, name="Gamman", value=1)
  setup_dg_equation(device, region_si, 'Lambda_e', 'Le', 'del_log_n1', 'del_log_n2', 'Le_eqn')
  setup_dg_contact(device, 'bot', 'Lambda_e', 'Le') 

def ActivateLh2():
  set_parameter(device=device, region=region_si, name="Gammap", value=3)
  set_parameter(device=device, region=region_ox, name="Gammap", value=1)
  setup_dg_equation(device, region_si, 'Lambda_h', 'Lh', 'del_log_p1', 'del_log_p2', 'Lh_eqn')
  setup_dg_contact(device, 'bot', 'Lambda_h', 'Lh') 

def set_default_models():
  for r in (region_si, region_ox):
    edge_from_node_model(device=device, region=r, node_model="x")
    CreateEdgeModel(device=device, region=r, model="xmid", expression="0.5*(x@n0 + x@n1)")
  set_parameter(name="T", value=300.0)
  CreateSolution(device=device, region=region_si, name="Potential")
  CreateSolution(device=device, region=region_ox, name="Potential")
  #
  # need Le and Lh for the Si Equation
  #
  setup_dg_variable(device, region_si, "Le")
  setup_dg_variable(device, region_si, "Lh")
  setup_dg_variable(device, region_ox, "Le")
  setup_dg_variable(device, region_ox, "Lh")

  CreateNodeModel(device=device, region=region_si, model="NetDoping", expression="Nconstant")

  #set the bias
  set_parameter(name=GetContactBiasName("bot"), value=0.0)
  set_parameter(name=GetContactBiasName("top"), value=0.0)
  #available solution reset
  node_solution(name='zero', device=device, region=region_ox)
  node_solution(name='zero', device=device, region=region_si)

def setup_potential_only():
  CreateOxidePotentialOnly(device=device, region=region_ox, update_type="default")
  CreateSiliconPotentialOnly(device=device, region=region_si)
  CreateSiliconOxideInterface(device=device, interface=interface_siox)
  CreateOxideContact(device=device, region=region_ox, contact=contact_ox)
  CreateSiliconPotentialOnlyContact(device=device, region=region_si, contact=contact_si)
  SetSiliconParameters(device=device, region=region_si)
  SetOxideParameters(device=device, region=region_ox)


device="MyDevice"
region_ox="MyOxRegion"
region_si="MySiRegion"
interface_siox="MySiOx"
contact_ox="top"
contact_si="bot"


try:
  tox = float(sys.argv[1])*1e-7
  carrier_var = sys.argv[2]
  pdfname = sys.argv[3]
except Exception as e:
  print '''args tox carrier pdfname
tox:     oxide thickness (nm)
carrier: Electrons or Holes
pdfname: pdf file name
'''
  raise



toxstring = r' t$_{ox}$=%s (nm)' % sys.argv[1]
xstring = r'distance from interface (nm)'

moscap.make_moscap(
  device=device,
  region_ox = region_ox,
  region_si = region_si,
  contact_ox = contact_ox,
  contact_si = contact_si,
  interface_siox = interface_siox,
  tox = tox,
  tsi = 60e-7,
  ox_material = 'Ox',
  si_material = 'Si'
)


#set_parameter(device=device, region=region_si, name="debug_level", value="verbose")
set_default_models()

setup_potential_only()

eqfuncs=None

def gamma_ramp_n():
  # ramp Gammanox
  ramp.rampparam(device=device, region=region_si, param="Gammanox", stop=1, step_size=0.1, min_step=0.001, solver_params=dg_solver_params)
  #ramp Gamman
  ramp.rampparam(device=device, region=region_si, param="Gamman",   stop=3.6, step_size=1, min_step=0.05, solver_params=dg_solver_params)

def gamma_ramp_p():
  # ramp Gammanox
  ramp.rampparam(device=device, region=region_si, param="Gammapox", stop=1, step_size=0.1, min_step=0.05, solver_params=dg_solver_params)
  #ramp Gamman
  ramp.rampparam(device=device, region=region_si, param="Gammap",   stop=3.6, step_size=0.1, min_step=0.05, solver_params=dg_solver_params)

if carrier_var == 'Electrons':
  dopings=['-1e17', '-1e18', '-1e19']
  dv=0.001
  factor=1e7 #from F/cm^2 to fF/um^2
  cstep=0.1
  qstep=0.1
  vstop=5
  vneg=-1.1
  eqfuncs=[ActivateLe2]
  gamma_ramps=[gamma_ramp_n]
elif carrier_var == 'Holes':
  #dopings=['1e19']
  dopings=['1e17', '1e18', '1e19']
  dv=0.001
  factor=1e7 #from F/cm^2 to fF/um^2
  cstep=0.1
  qstep=0.1
  vstop=-5
  vneg=1.1
  eqfuncs=[ActivateLh2]
  gamma_ramps=[gamma_ramp_p]


classical_data = [None]*len(dopings)

solver_params = {
  'type' : "dc",
  'relative_error' : 1e-11,
  'absolute_error' : 1,
  'maximum_iterations' : 50
}
for index, d in enumerate(dopings):
  cdict = {
    'doping' : format_doping(float(d)),
    'q' : [],
    'v' : [],
  }


  set_parameter(device=device, region=region_si, name="Nconstant", value=float(d))
  set_node_values(device=device, region=region_si, name="Potential", init_from="zero")
  set_node_values(device=device, region=region_ox, name="Potential", init_from="zero")

  set_parameter(name=GetContactBiasName("top"), value=vneg)
  devsim.solve(**solver_params)

  def mycb():
    res = simulate_charge(device=device, contact="top", equation="PotentialEquation", solver_params=solver_params)
    cdict['v'].append(res[0])
    cdict['q'].append(res[1])
    
  ramp.rampparam(device="", region="", param=GetContactBiasName("top"), stop=vstop, step_size=cstep, min_step=0.01, solver_params=solver_params, callback=mycb)

  cdict['Electrons'] = numpy.array(get_node_model_values(device=device, region=region_si, name="IntrinsicElectrons"))
  cdict['Holes']     = numpy.array(get_node_model_values(device=device, region=region_si, name="IntrinsicHoles"))
  cdict['Potential'] = numpy.array(get_node_model_values(device=device, region=region_si, name="Potential"))
  cdict['q'] = numpy.array(cdict['q'])
  cdict['v'] = numpy.array(cdict['v'])
  classical_data[index] = cdict


# The new dg physics
#
set_dg_parameters(device, region_si)
set_dg_parameters(device, region_ox)

setup_atox(device, region_si)
setup_dg_si_potential_only(device, region_si)
setup_dg_ox(device, region_ox)

for eq in eqfuncs:
  eq()

quantum_data = [None]*len(dopings)

dg_solver_params = {
  'type' : "dc",
  'relative_error' : 1e-5,
  'absolute_error' : 1,
  'maximum_iterations' : 50
}


for index, d in enumerate(dopings):
  qdict = {
    'doping' : format_doping(float(d)),
    'q' : [],
    'v' : [],
  }

  set_parameter(device=device, region=region_si, name="Nconstant", value=float(d))
  set_node_values(device=device, region=region_si, name="Potential", init_from="zero")
  set_node_values(device=device, region=region_si, name="Le", init_from="zero")
  set_node_values(device=device, region=region_si, name="Lh", init_from="zero")
  set_node_values(device=device, region=region_ox, name="Potential", init_from="zero")

  set_parameter(name=GetContactBiasName("top"), value=0.0)
  print d
  set_parameter(device=device, region=region_si, name="Gamman", value=10)
  set_parameter(device=device, region=region_si, name="Gammanox", value=1)
  set_parameter(device=device, region=region_si, name="Gammap", value=0.1)
  set_parameter(device=device, region=region_si, name="Gammapox", value=0.00001)
  devsim.solve(**dg_solver_params)

  for gr in gamma_ramps:
    gr()


  # ramp top bias negative
  ramp.rampparam(device="", region="", param=GetContactBiasName("top"), stop=vneg, step_size=0.1, min_step=0.01, solver_params=dg_solver_params, callback=None)


  def mycb():
    res = simulate_charge(device=device, contact="top", equation="PotentialEquation", solver_params=dg_solver_params)
    qdict['v'].append(res[0])
    qdict['q'].append(res[1])
    
  ramp.rampparam(device="", region="", param=GetContactBiasName("top"), stop=vstop, step_size=qstep, min_step=0.01, solver_params=dg_solver_params, callback=mycb)

  qdict['Electrons'] = numpy.array(get_node_model_values(device=device, region=region_si, name="IntrinsicElectrons"))
  qdict['Holes']     = numpy.array(get_node_model_values(device=device, region=region_si, name="IntrinsicHoles"))
  qdict['Potential'] = numpy.array(get_node_model_values(device=device, region=region_si, name="Potential"))
  qdict['Le']        = numpy.array(get_node_model_values(device=device, region=region_si, name="Le"))
  qdict['Lh']        = numpy.array(get_node_model_values(device=device, region=region_si, name="Lh"))
  qdict['q'] = numpy.array(qdict['q'])
  qdict['v'] = numpy.array(qdict['v'])
  quantum_data[index] = qdict

#print max(get_node_model_values(device=device, region=region_si, name='n_classical'))
#print max(get_node_model_values(device=device, region=region_si, name='n_quantum'))

if carrier_var == 'Electrons':
  edge_model(device=device, region=region_si, name="surface_term", equation="del_log_n1*EdgeCouple")
  edge_model(device=device, region=region_si, name="edge_volume_term", equation="del_log_n2*EdgeNodeVolume")
  node_model(device=device, region=region_si, name="volume_term", equation="Le_eqn*NodeVolume")
elif carrier_var == 'Holes':
  edge_model(device=device, region=region_si, name="surface_term", equation="del_log_p1*EdgeCouple")
  edge_model(device=device, region=region_si, name="edge_volume_term", equation="del_log_p2*EdgeNodeVolume")
  node_model(device=device, region=region_si, name="volume_term", equation="Lh_eqn*NodeVolume")
else:
  raise RuntimeError("Unknown Carrier Type")

figs = []
fig = debug_plot(device=device, region=region_si)
fig.suptitle(('debugging plot for %s,' + toxstring) % (qdict['doping'],) )
figs.append(fig)


##
##plt.plot(
##  get_edge_model_values(device=device, region=region_si, name='xmid'),
##  get_edge_model_values(device=device, region=region_si, name='surface_term'), 'x',
##  label="surf"
##)
###plt.show()
##plt.plot(
##  get_edge_model_values(device=device, region=region_si, name='xmid'),
##  get_edge_model_values(device=device, region=region_si, name='edge_volume_term'), '+',
##  label="vol edge"
##)
##plt.plot(
##  get_node_model_values(device=device, region=region_si, name='x'),
##  get_node_model_values(device=device, region=region_si, name='volume_term'), 'o',
##  label="vol node"
##)
##plt.legend(loc='best')
##plt.show()
##
##
###plt.ylim(0,4e20)
###plt.show()
##plt.plot(
##  get_node_model_values(device=device, region=region_si, name='x'),
##  get_node_model_values(device=device, region=region_si, name='Le'), '+',
##  label="classical"
##)
###plt.plot(
###  get_node_model_values(device=device, region=region_si, name='x'),
###  get_node_model_values(device=device, region=region_si, name='potential_effective'),
###  label="effective potential"
###)
###plt.xlim(0, 7e-7)
###plt.ylim(1e18,1e22)
##plt.show()
##
#x=numpy.array(get_node_model_values(device=device, region=region_si, name='x'))
#for i in range(len(dopings)):
#  plt.plot(
#    x,
#    classical_data[i]['Electrons'],
#    label="Classical"
#  )
#  plt.plot(
#    x,
#    quantum_data[i]['Electrons'],python_packages.
#    label="Density Gradient"
#  )
#  plt.title(classical_data[i]['doping'])
#  plt.xlim(0, 7e-7)
#  plt.ylim(0,1e21)
#  plt.legend(loc='best')
#  plt.show()
#
x=abs(numpy.array(get_node_model_values(device=device, region=region_si, name='x')))*1e7
for i in range(len(dopings)):
  fig = plt.figure()
  plt.plot(
    x,
    classical_data[i][carrier_var],
    label="Classical"
  )
  plt.plot(
    x,
    quantum_data[i][carrier_var],
    label="Density Gradient"
  )
  plt.title(('%s' + toxstring) % (classical_data[i]['doping'],))
  plt.xlim(0, 7)
  plt.ylim(0,1e21)
  plt.legend(loc='best')
  plt.xlabel(xstring)
  plt.ylabel('%s (#/cm$^3$)' % (carrier_var,))
  #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False, useMathText=True)
  ax = plt.gca()
  ax.yaxis.set_major_formatter(y_formatter)
  figs.append(fig)
#
#
#
colors=('k', 'b', 'g')
fig = plt.figure()
for i in range(len(dopings)):
  plt.plot(classical_data[i]['v'], classical_data[i]['q'], '%s--' % colors[i], label='%s' % classical_data[i]['doping'])
  plt.plot(quantum_data[i]['v'],   quantum_data[i]['q'],   '%s-' % colors[i], label='DG')
plt.title('CV Curves' + toxstring)
plt.xlabel(r'V$_g$ (V)')
plt.ylabel(r'C (fF/$\mu$m)')
plt.legend(loc='best')
figs.append(fig)

with PdfPages(pdfname) as pdf:
  for f in figs:
    pdf.savefig(f)
    plt.close(f)

