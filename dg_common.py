from dg_physics import *

def set_dg_parameters(device, region):
    ##From Cogenda default Si physics:
    #TODO: parameter units
    set_parameter(device=device, region=region, name="h", value=6.62606876e-34)  # J*s
    set_parameter(device=device, region=region, name="hbar", value=1.054571596e-34) # J*s
    set_parameter(device=device, region=region, name="Gamman", value=1)
    set_parameter(device=device, region=region, name="Gammap", value=1)
    set_parameter(device=device, region=region, name="Gammanox", value=1)
    set_parameter(device=device, region=region, name="Gammapox", value=1)
    set_parameter(device=device, region=region, name="mass_electron", value=9.10938188e-31) # kg
    node_model(device=device, region=region, name="meff_elec", equation="1.0903*mass_electron")
    node_model(device=device, region=region, name="meff_hole", equation="1.1525*mass_electron")
    node_model(device=device, region=region, name="b_n", equation="(1e4*Gamman*hbar^2)/(q*6*meff_elec)") # eV*cm^2
    node_model(device=device, region=region, name="b_p", equation="(1e4*Gammap*hbar^2)/(q*6*meff_hole)") # eV*cm^2
    # Oxide
    #
    #node_model(device=device, region=region, name="b_nox", equation="1e-100") # eV*cm^2
    node_model(device=device, region=region, name="b_nox", equation="1e4*Gammanox*hbar^2/(6*q*0.14*mass_electron)") # eV*cm^2
    node_model(device=device, region=region, name="b_pox", equation="1e4*Gammapox*hbar^2/(6*q*1.0*mass_electron)") # eV*cm^2
    #node_model(device=device, region=region, name="b_nox", equation="10*hbar^2/(6*q*0.14*mass_electron)") # eV*cm^2
    #node_model(device=device, region=region, name="b_pox", equation="10*hbar^2/(6*q*0.14*mass_electron)") # eV*cm^2

    #node_model(device=device, region=region, name="b_nox", equation="(1e4*Gammanox*hbar^2)/(q*6*meff_elec)") # eV*cm^2
    #node_model(device=device, region=region, name="b_pox", equation="(1e4*Gammapox*hbar^2)/(q*6*meff_hole)") # eV*cm^2
    node_model(device=device, region=region, name="x_np", equation="1e2*hbar/(2*0.4*mass_electron*q*3.15)^0.5") # cm
    node_model(device=device, region=region, name="x_pp", equation="1e2*hbar/(2*0.4*mass_electron*q*4.50)^0.5") # cm
    #const PetscScalar b_nox = hbar*hbar/(6*e*0.14*me);
    #const PetscScalar b_pox = hbar*hbar/(6*e*1.0*me);
    #const PetscScalar bn = mt->band->Gamman()*hbar*hbar/(6*e*mt->band->EffecElecMass(T));
    #const PetscScalar bp = mt->band->Gammap()*hbar*hbar/(6*e*mt->band->EffecHoleMass(T));
    #  const PetscScalar x_np = hbar/sqrt(2*0.4*me*3.15*eV);
    #  const PetscScalar x_pp = hbar/sqrt(2*0.4*me*4.50*eV);
    #    eV = 1.602176462e-19*J;
    #    e    = 1.602176462e-19*C;
    #    ELECMASS  = 1.0903*me;
    #    HOLEMASS  = 1.1525*me;

    #
    # Edge Integration
    #
    em = (
        ("b_n", "b_ne", "arithmetic", ()),
      ("b_p", "b_pe", "arithmetic", ()),
    )
    for e in em:
        edge_average_model(device=device, region=region, node_model=e[0], edge_model=e[1], average_type=e[2])
        for d in e[3]:
            edge_average_model(device=device, region=region, node_model=e[0], edge_model=e[1], average_type=e[2], derivative=d)


def setup_atox(device, region):
    CreateSolution(device=device, region=region, name="AtOx")
    for i in get_interface_list(device=device):
        rlist = get_region_list(device=device, interface=i)  
        if region in rlist:
            if rlist[0] == region:
                interface_model(device=device, interface=i, name="iindex", equation="node_index@r0")
            else:
                interface_model(device=device, interface=i, name="iindex", equation="node_index@r1")
            nlist = get_interface_model_values(device=device, interface=i, name="iindex")
            for j in nlist:
                print(j)
                set_node_value(device=device, region=region, index=int(j), name="AtOx", value=1.0)

def setup_dg_variable(device, region, variable):
    CreateSolution(device, region, variable)
    edge_from_node_model(device=device, region=region, node_model=variable)

def setup_log_si(device, region):
    #
    # log_n, log_p
    #
    em = (
        ("d_l_n",           "log(Electrons@n0/Electrons@n1)*EdgeInverseLength"),
      ("d_l_n:Electrons@n0", "(Electrons@n0^(-1))*EdgeInverseLength"),
      ("d_l_n:Electrons@n1", "-(Electrons@n1^(-1))*EdgeInverseLength"),
      ("d_l_p",           "log(Holes@n0/Holes@n1)*EdgeInverseLength"),
      ("d_l_p:Holes@n0", "(Holes@n0^(-1))*EdgeInverseLength"),
      ("d_l_p:Holes@n1", "-(Holes@n1^(-1))*EdgeInverseLength"),
      # This is since we need to wipe out previous derivative
      ("d_l_n:Potential@n0", "0"),
      ("d_l_n:Potential@n1", "0"),
      ("d_l_n:Le@n0", "0"),
      ("d_l_n:Le@n1", "0"),
      ("d_l_p:Potential@n0", "0"),
      ("d_l_p:Potential@n1", "0"),
      ("d_l_p:Lh@n0", "0"),
      ("d_l_p:Lh@n1", "0"),
    )
    for e in em:
        edge_model(device=device, region=region, name=e[0], equation=e[1])

def setup_log_ox(device, region):

    edge_from_node_model(device=device, region=region, node_model="Le")
    edge_from_node_model(device=device, region=region, node_model="Lh")

    #
    # Need the V_t models
    #
    CreateVT(device, region, ())

    #
    # log_n, log_p
    #
    em = (
        # could use gradient model
      ("d_l_n",           "EdgeInverseLength*(Le@n1 - Le@n0)/V_t_edge"),
      ("d_l_n:Le@n0",     "-EdgeInverseLength/V_t_edge"),
      ("d_l_n:Le@n1",     "EdgeInverseLength/V_t_edge"),
      #("d_l_n:Potential@n0", "0"),
      #("d_l_n:Potential@n1", "0"),
      #("d_l_n:Le@n0", "0"),
      #("d_l_n:Le@n1", "0"),

      ("d_l_p",           "EdgeInverseLength*(Lh@n1 - Lh@n0)/V_t_edge"),
      ("d_l_p:Lh@n0",     "-EdgeInverseLength/V_t_edge"),
      ("d_l_p:Lh@n1",     "EdgeInverseLength/V_t_edge"),
      #("d_l_p:Potential@n0", "0"),
      #("d_l_p:Potential@n1", "0"),
      #("d_l_p:Lh@n0", "0"),
      #("d_l_p:Lh@n1", "0"),
    )
    for e in em:
        edge_model(device=device, region=region, name=e[0], equation=e[1])

def setup_log_si_potential_only(device, region):

    edge_from_node_model(device=device, region=region, node_model="Le")
    edge_from_node_model(device=device, region=region, node_model="Lh")
    edge_from_node_model(device=device, region=region, node_model="NIE")

    #
    # Need the V_t models
    #
    CreateVT(device, region, ())

    #
    # log_n, log_p
    #
    em = (
        # could use gradient model
      #("d_l_n",           "EdgeInverseLength*((Le@n1 - Le@n0)/V_t_edge)"),
      #("d_l_n",           "EdgeInverseLength*((-V_t_edge*log(NIE_edge) + (Le@n1 - Le@n0))/V_t_edge)"),
      #("d_l_n",           "EdgeInverseLength*((-V_t_edge*log(NIE_edge) + (Potential@n0 - Potential@n1))/V_t_edge)"),
      #TODO: add band edge element later (NIE)
      ("d_l_n",           "EdgeInverseLength*(((Potential@n1 - Potential@n0) + (Le@n0 - Le@n1))/V_t_edge)"),
      #("d_l_n",           "EdgeInverseLength*((V_t_edge*log(NIE@n0/NIE@n1) + (Potential@n0 - Potential@n1) + (Le@n1 - Le@n0))/V_t_edge)"),
      #("d_l_n",           "EdgeInverseLength*((V_t_edge*log(NIE_edge) + (Potential@n0 - Potential@n1) + (Le@n1 - Le@n0))/V_t_edge)"),
      ("d_l_n:Potential@n0",     "-EdgeInverseLength/V_t_edge"),
      ("d_l_n:Potential@n1",     "EdgeInverseLength/V_t_edge"),
      ("d_l_n:Le@n0",            "EdgeInverseLength/V_t_edge"),
      ("d_l_n:Le@n1",            "-EdgeInverseLength/V_t_edge"),
      #("d_l_p",           "EdgeInverseLength*((Lh@n1 - Lh@n0)/V_t_edge)"),
      ("d_l_p",           "EdgeInverseLength*(((Potential@n0 - Potential@n1) + (Lh@n0 - Lh@n1))/V_t_edge)"),
      #("d_l_p",           "EdgeInverseLength*((Potential@n1 - Potential@n0) + (Lh@n1 - Lh@n0))/V_t_edge"),
      ("d_l_p:Potential@n0",     "EdgeInverseLength/V_t_edge"),
      ("d_l_p:Potential@n1",     "-EdgeInverseLength/V_t_edge"),
      ("d_l_p:Lh@n0",     "EdgeInverseLength/V_t_edge"),
      ("d_l_p:Lh@n1",     "-EdgeInverseLength/V_t_edge"),
    )
    for e in em:
        edge_model(device=device, region=region, name=e[0], equation=e[1])

def setup_dg_edge(device, region, carrier_name, parameter_name, variables):
    #
    # Looking at the source code, it appears that flux is added into node 0 and subtracted from node 1
    # and volume edge terms are added to node 0 and added to node 1
    # in the future we may change the source to reverse the sign on the volume edge term
    #
    #
    # Volume integration term
    #
    model_name="del_log_%s1" % (carrier_name,)
    expression="0.5*%s*d_l_%s" % (parameter_name, carrier_name)

    edge_model(device=device, region=region, name=model_name, equation=expression)
    for v in variables:
        m0 = ":%s@n0" % (v,)
        m1 = ":%s@n1" % (v,)
        edge_model(device=device, region=region, name=(model_name + m0), equation=(expression + m0))
        edge_model(device=device, region=region, name=(model_name + m1), equation=(expression + m1))

    model_name="del_log_%s2" % (carrier_name,)
    expression="0.25*%s*1e8*(1e-4*d_l_%s)^2" % (parameter_name, carrier_name)
    edge_model(device=device, region=region, name=model_name, equation=expression)
    expression="0.50*%s*d_l_%s*d_l_%s" % (parameter_name, carrier_name, carrier_name)
    for v in variables:
        m0 = ":%s@n0" % (v,)
        m1 = ":%s@n1" % (v,)
        edge_model(device=device, region=region, name=(model_name + m0), equation=(expression + m0))
        edge_model(device=device, region=region, name=(model_name + m1), equation=(expression + m1))

def setup_dg_volume(device, region, normalized=False):
    #
    # Lambda Equations
    #
    #em = (
    #  ("Le_eqn",         "Le/b_n"),
    #  ("Le_eqn:Le",      "1/b_n"),
    #  ("Lh_eqn",         "Lh/b_p"),
    #  ("Lh_eqn:Lh",      "1/b_p"),
    #)
    em = (
        #    ("Le_eqn",         "Le - AtOx*1e2"),
        ("Le_eqn",         "Le - AtOx*b_nox*SurfaceArea/(NodeVolume*x_np)"),
      ("Le_eqn:Le",      "1"),
        #    ("Lh_eqn",         "Lh"),
      ("Lh_eqn",         "Lh - AtOx*b_pox*SurfaceArea/(NodeVolume*x_pp)"),
      ("Lh_eqn:Lh",      "1"),
    )

    for e in em:
        node_model(device=device, region=region, name=e[0], equation=e[1])

#
# TODO: homojunction
#

def setup_edge_volume(device, region):
    #
    # Set edge volumes for integration
    #
    if get_dimension(device=device) == 1:
        edge_model(device=device, region=region, name="EdgeNodeVolume", equation="0.5*EdgeCouple*EdgeLength")
    elif get_dimension(device=device) == 2:
        edge_model(device=device, region=region, name="EdgeNodeVolume", equation="0.25*EdgeCouple*EdgeLength")
    else:
        raise NameError("Unhandled dimension")
    set_parameter(device=device, region=region, name="edge_node0_volume_model", value="EdgeNodeVolume")
    set_parameter(device=device, region=region, name="edge_node1_volume_model", value="EdgeNodeVolume")


def setup_dg_si(device, region):
    setup_edge_volume(device, region)

    setup_log_si(device, region)

    setup_dg_edge(device, region, 'n', 'b_ne', ('Electrons',))
    setup_dg_edge(device, region, 'p', 'b_pe', ('Holes',))

    setup_dg_volume(device, region)


def setup_dg_si_potential_only(device, region):
    setup_edge_volume(device, region)

    setup_log_si_potential_only(device, region)

    setup_dg_edge(device, region, 'n', 'b_ne', ('Potential', 'Le'))
    setup_dg_edge(device, region, 'p', 'b_pe', ('Potential', 'Lh'))

    setup_dg_volume(device, region)


def setup_dg_ox(device, region):
    setup_edge_volume(device, region)

    setup_log_ox(device, region)

    setup_dg_edge(device, region, 'n', 'b_ne', ('Le',))
    setup_dg_edge(device, region, 'p', 'b_pe', ('Lh',))

    setup_dg_volume(device, region)

def setup_dg_equation(device, region, e_name, v_name, e_m, e_v_m, n_m):
    equation(device=device, region=region, name=e_name, variable_name=v_name,
             edge_model=e_m, edge_volume_model=e_v_m, node_model=n_m, variable_update="default")

# TODO: generalize multiple contacts by testing existence
def setup_dg_contact(device, contact, eq_name, v_name):
    #TODO: test if these models already exist on the region
    m = eq_name + '_' + contact
    CreateContactNodeModel(device, contact, m, v_name)
    d = m + ':' + v_name
    CreateContactNodeModel(device, contact, d, '1')

    contact_equation(device=device , contact=contact, name=eq_name, variable_name=v_name,
                     node_model=m)

def setup_dg_interface(device, interface, eq_name, v_name):
    model_name = CreateContinuousInterfaceModel(device, interface, v_name)
    interface_equation(device=device, interface=interface, name=eq_name, variable_name=v_name, interface_model=model_name, type="continuous")

