# Make cos(theta*) and w_a,w_s from the reco t and tbar p4
from ROOT import *
from math import sqrt
import math

# Get cos(theta*) for reco top pair p4
def get_angles(reco_t,reco_tbar):
    
    #Lorentz Declartions

    sqrt_s=8000;
    beam_energy=sqrt_s/2;

    S = TLorentzRotation();

    #intialize the rotations to the identity
    R_data = TLorentzRotation() ;

    # Set up t and tbar p4
    Top_Data = reco_t.Clone()
    ATop_Data = reco_tbar.Clone()  
    
    #Make the 4-vector of the ttbar
    Q_Data = Top_Data + ATop_Data;
    ttbar_mass_data=Q_Data.Mag(); # to return Mtt

    #defining the Px, Py,and Pz, and energies to boost into the ttbar rest frame
    Bx_data = -1*Q_Data.Px()/Q_Data.E();  
    By_data = -1*Q_Data.Py()/Q_Data.E();  
    Bz_data = -1*Q_Data.Pz()/Q_Data.E();
    Qt = sqrt(Q_Data.Px()*Q_Data.Px()+Q_Data.Py()*Q_Data.Py());

    #Feynman x
    x_f = 2*Q_Data.Pz()/sqrt_s; # to return x_f
    
    #Creating the Lorentz 4 vectors of the protons
    Proton_data = TLorentzVector(0,0,sqrt(beam_energy*beam_energy -1*1),beam_energy)
    Proton2_data = TLorentzVector(0,0,-1*Proton_data.Pz(),beam_energy)
    
    #Doing the boost
    R_data = R_data.Boost(Bx_data,By_data,Bz_data);
    Top_Data = R_data*Top_Data;
    ATop_Data = R_data*ATop_Data;
    Proton_data = R_data*Proton_data;
    Proton2_data = R_data*Proton2_data;
    #Reset the boost
    R_data=S;
    #Define three vectors for P,Pbar,top, quark and antiquark in ttbar c.m frame
    top_vec = Top_Data.Vect()
    proton_vec = Proton_data.Vect()
    proton2_vec = Proton2_data.Vect()
    
    #Flip the larger one between proton and proton2, and flip the antiquark direction
    if proton_vec.Mag()>proton2_vec.Mag() : proton_vec=-1.0*proton_vec;
    else : proton2_vec=-1.0*proton2_vec;
    
    #Normalize vectors
    top_vec = top_vec*(1.0/top_vec.Mag());
    proton_vec = proton_vec*(1.0/proton_vec.Mag());
    proton2_vec = proton2_vec*(1.0/proton2_vec.Mag());
    #find the unit bisectors
    bisector_vec = (proton_vec+proton2_vec)*(1.0/(proton_vec+proton2_vec).Mag());
    #find the CS angle
    cos_theta_cs=top_vec.Dot(bisector_vec);

    return [x_f,ttbar_mass_data,cos_theta_cs]



<<<<<<< HEAD
# Get cos(theta*) for reco top pair p4
def get_angles_v1(reco_t,reco_tbar):
    
    #Lorentz Declartions

    sqrt_s=8000;
    beam_energy=sqrt_s/2;

    S = TLorentzRotation();

    #intialize the rotations to the identity
    R_data = TLorentzRotation() ;
      
    #Vectors of top and antitop
    Top_Data = TLorentzVector();
    ATop_Data = TLorentzVector() ;

    # Set up t and tbar p4
    Top_Data = reco_t.Clone()
    ATop_Data = reco_tbar.Clone()

    #assign the quark and antiquark fourvectors from the last two particles in the list of mc_ particles
    quark = TLorentzVector(0.0,0.0,sqrt(beam_energy*beam_energy -1*1),beam_energy)
    antiquark = TLorentzVector(0.0,0.0,-1*quark.Pz(),beam_energy)   
    
    #Make the 4-vector of the ttbar
    Q_Data = Top_Data + ATop_Data;
    ttbar_mass_data=Q_Data.Mag(); # to return Mtt

    #defining the Px, Py,and Pz, and energies to boost into the ttbar rest frame
    Bx_data = -1*Q_Data.Px()/Q_Data.E();  
    By_data = -1*Q_Data.Py()/Q_Data.E();  
    Bz_data = -1*Q_Data.Pz()/Q_Data.E();
    Qt = sqrt(Q_Data.Px()*Q_Data.Px()+Q_Data.Py()*Q_Data.Py());
    beta_mc = (Bx_data,By_data,Bz_data)

    #Feynman x
    x_f = 2*Q_Data.Pz()/sqrt_s; # to return x_f
    
    #Creating the Lorentz 4 vectors of the protons
    Proton_data = TLorentzVector(0,0,sqrt(beam_energy*beam_energy -1*1),beam_energy)
    Proton2_data = TLorentzVector(0,0,-1*Proton_data.Pz(),beam_energy)
    
    #Doing the boost
    R_data = R_data.Boost(Bx_data,By_data,Bz_data);
    Top_Data = R_data*Top_Data;
    ATop_Data = R_data*ATop_Data;
    Proton_data = R_data*Proton_data;
    Proton2_data = R_data*Proton2_data;
    proton1_e = Proton_data.E()
    proton2_e = Proton2_data.E()

    #Reset the boost
    R_data=S.Clone();
    #Define three vectors for P,Pbar,top, quark and antiquark in ttbar c.m frame
    top_data = TVector3(Top_Data.Px(),Top_Data.Py(),Top_Data.Pz())
    proton_data = TVector3(Proton_data.Px(),Proton_data.Py(),Proton_data.Pz())
    proton2_data = TVector3(Proton2_data.Px(),Proton2_data.Py(),Proton2_data.Pz())
   
    # study the consequenc of boost of proton momentum
    delta_mag = proton_data.Mag()
    beta1 = proton_data.Clone()
    beta1 *= 1/proton1_e
    beta1 = (beta1.Px(),beta1.Py(),beta1.Pz())
    beta2 = proton2_data.Clone()
    beta2 *= 1/proton2_e
    beta2 = (beta2.Px(),beta2.Py(),beta2.Pz())


    #Flip the larger one between proton and proton2, and flip the antiquark direction    
    if proton_data.Mag()>proton2_data.Mag() : proton_data=-1.0*proton_data;
    else : proton2_data=-1.0*proton2_data;
    
    #Normalize vectors
    top_data = top_data*(1.0/top_data.Mag());
    proton_data = proton_data*(1.0/proton_data.Mag());
    proton2_data = proton2_data*(1.0/proton2_data.Mag());
    #find the unit bisectors
    bisector_data = (proton_data+proton2_data)*(1.0/(proton_data+proton2_data).Mag());
    #find the CS angle
    cos_theta_cs_data=top_data.Dot(bisector_data);
    #find the angle bewtween unit bisector and z axis
    vec_z = TVector3(0,0,1)
    cs_bisector = vec_z.Dot(bisector_data)

    return [x_f,ttbar_mass_data,cos_theta_cs_data,cs_bisector,beta_mc,beta1,beta2,delta_mag]




# Get cos(theta*) for top pair p4 using mc p4 and mc q and qbar
def get_true_angles(reco_t,reco_tbar,q_p4,qbar_p4):
=======
# Get true cos(theta*) using mc tops p4 and mc q z direction. If no qz direction provided, say, in gg,q1q2 case, use regular cs direction def.
def get_true_angles(mc_t,mc_tbar,qz=0):
>>>>>>> mydev
    
    #Lorentz Declartions

    sqrt_s=8000;
    beam_energy=sqrt_s/2;

    S = TLorentzRotation();

    #intialize the rotations to the identity
    R_data = TLorentzRotation() ;

    # Set up t and tbar p4
    Top_Data = mc_t.Clone()
    ATop_Data = mc_tbar.Clone()  
    
    #Make the 4-vector of the ttbar
    Q_Data = Top_Data + ATop_Data;
    ttbar_mass_data=Q_Data.Mag(); # to return Mtt

    #defining the Px, Py,and Pz, and energies to boost into the ttbar rest frame
    Bx_data = -1*Q_Data.Px()/Q_Data.E();  
    By_data = -1*Q_Data.Py()/Q_Data.E();  
    Bz_data = -1*Q_Data.Pz()/Q_Data.E();
    Qt = sqrt(Q_Data.Px()*Q_Data.Px()+Q_Data.Py()*Q_Data.Py());

    #Feynman x
    x_f = 2*Q_Data.Pz()/sqrt_s; # to return x_f
    
    #Creating the Lorentz 4 vectors of the protons
    Proton_data = TLorentzVector(0,0,sqrt(beam_energy*beam_energy -1*1),beam_energy)
    Proton2_data = TLorentzVector(0,0,-1*Proton_data.Pz(),beam_energy)
    
    #Doing the boost
    R_data = R_data.Boost(Bx_data,By_data,Bz_data);
    Top_Data = R_data*Top_Data;
    ATop_Data = R_data*ATop_Data;
    Proton_data = R_data*Proton_data;
    Proton2_data = R_data*Proton2_data;
    #Reset the boost
    R_data=S;
    #Define three vectors for P,Pbar,top, quark and antiquark in ttbar c.m frame
    top_vec = Top_Data.Vect()
    proton_vec = Proton_data.Vect()
    proton2_vec = Proton2_data.Vect()
    
    #Flip the larger one between proton and proton2, and flip the antiquark direction
    if proton_vec.Mag()>proton2_vec.Mag() : proton_vec=-1.0*proton_vec;
    else : proton2_vec=-1.0*proton2_vec;
    
    #Normalize vectors
    top_vec = top_vec*(1.0/top_vec.Mag());
    proton_vec = proton_vec*(1.0/proton_vec.Mag());
    proton2_vec = proton2_vec*(1.0/proton2_vec.Mag());
    # find the unit bisectors direction
    bisector_vec = (proton_vec+proton2_vec)*(1.0/(proton_vec+proton2_vec).Mag());
    # Make sure unit bisectors direction is the same as mc_q direction
    if bisector_vec*qz < -0.1 : # aka bisector and true qz are in opposite direction and qz exists!
        bisector_vec *= -1.0
    #find the CS angle
<<<<<<< HEAD
    cos_theta_cs_true=top_data.Dot(bisector_mc);

    return [x_f,ttbar_mass_data,cos_theta_cs_true]

# Get cos(theta*) for top pair p4 using mc p4 and mc init particles
# Use quark direction as positive direction
def get_true_angles_v1(reco_t,reco_tbar,q_pz,qbar_pz,q_id,qbar_id):
    
    #Lorentz Declartions

    sqrt_s=8000;
    beam_energy=sqrt_s/2;

    S = TLorentzRotation();

    #intialize the rotations to the identity
    R_data = TLorentzRotation() ;
      
    #Vectors of beam1 and beam2 for unit bisector purpose
    # TLorentzVector(px,py,pz,E)
    beam1 = TLorentzVector(0.0,0.0,sqrt(beam_energy*beam_energy -1*1),beam_energy)
    beam2 = TLorentzVector(0.0,0.0,-1*beam1.Pz(),beam_energy)

    # Set up t and tbar p4
    Top_Data = reco_t.Clone()
    ATop_Data = reco_tbar.Clone()

    # check the cos_theta in lab frame
    top_vec = Top_Data.Vect()
    top_vec *= 1/top_vec.Mag()
    z_vec = TVector3(0,0,1)
    cs_lab = top_vec.Dot(z_vec)

    # find the true quark direction
    q_pz_list = [q_pz,qbar_pz]
    q_id_list = [q_id,qbar_id]
    positive_z = 0
    num_qs = 0
    for i in range(2):
        if q_id_list[i] in [1,2,3,4]:
            num_qs += 1
            if q_pz_list[i]>0:
                positive_z = 1
            else:
                positive_z = -1
    if num_qs != 1 :
        positive_z = 0

  
    #Make the 4-vector of the ttbar
    Q_Data = Top_Data + ATop_Data;
    ttbar_mass_data=Q_Data.Mag(); # to return Mtt
=======
    cos_theta_cs=top_vec.Dot(bisector_vec);
>>>>>>> mydev

    return [x_f,ttbar_mass_data,cos_theta_cs]

<<<<<<< HEAD
    #Feynman x
    # x_f = (quark.E()-antiquark.E())/beam_energy
    x_f = 2*Q_Data.Pz()/sqrt_s; # to return x_f 

    #Doing the boost
    R_data = R_data.Boost(Bx_data,By_data,Bz_data);
    Top_Data = R_data*Top_Data;
    ATop_Data = R_data*ATop_Data;
    beam1 = R_data*beam1
    beam2 = R_data*beam2
    #Reset the boost
    R_data=S.Clone();
    #Define three vectors for P,Pbar,top, quark and antiquark in ttbar c.m frame
    top_data = TVector3(Top_Data.Px(),Top_Data.Py(),Top_Data.Pz())
    top_data = top_data*(1.0/top_data.Mag());

    beam1_vec = beam1.Vect()
    beam2_vec = beam2.Vect()
    if beam1_vec.Mag()>beam2_vec.Mag():
        beam1_vec *= -1.0
    else:
        beam2_vec *= -1.0
    bisector = beam1_vec*(1.0/beam1_vec.Mag())+beam2_vec*(1.0/beam2_vec.Mag())

    # Find true quark direction
    if positive_z*bisector.Pz()< -0.1:
        bisector *= -1.0
    true_quark_direction = bisector*(1.0/bisector.Mag())
    #find the CS angle
    cos_theta_cs_true=top_data.Dot(true_quark_direction);

    return [cos_theta_cs_true,cs_lab]
=======
def Get_beta(m1,m2,mtt):
    return math.sqrt(1-2*(m1*m1+m2*m2)/pow(mtt,2)+pow((m1*m1-m2*m2),2)/pow(mtt,4))

>>>>>>> mydev

# Get reweighting weights for qqbar->ttbar templates
def GetAnglesWeights(Top_MC,ATop_MC,cos_theta_cs_mc,alpha_input=-0.129):
    m1   = Top_MC.M();
    m2   = ATop_MC.M();
    mtt = (Top_MC+ATop_MC).M()
    beta_mc   = Get_beta(m1,m2,mtt)
    if TMath.IsNaN(beta_mc) : 
        print 'Warning! beta_mc is negative!'
        return ()
    # weights (cf equation 5 in note)
    alpha = 0.0;
    #      THESE ALPHA VALUES MUST BE CHANGED FOR GENERATORS OTHER THAN MADGRAPH 5
    # if (nValidJets[0] == 4)
    #   alpha = -0.228;
    # else if (nValidJets[0] == 5)
    #   alpha = 0.010;
    # VALUES FOR CT10 POWHEG
    # if (nValidJets == 4)
    #   alpha = -0.256;
    # else if (nValidJets == 5)
    #   alpha = 0.143;
    # combined average alpha
    alpha = alpha_input;
    one_m_b2 = 1.0-beta_mc*beta_mc;
    b2c2 = beta_mc*beta_mc*cos_theta_cs_mc*cos_theta_cs_mc;
    otb2 = (1.0/3.0)*beta_mc*beta_mc;
    denom = 1.0+b2c2+one_m_b2+alpha*(1.0-b2c2);
    w_a = 2.0 * ((1.0+otb2+one_m_b2+alpha*(1.0-otb2))/denom) * cos_theta_cs_mc; 
    w_s_xi = one_m_b2/denom;
    w_a_xi = 2.0*(one_m_b2/denom)*cos_theta_cs_mc;
    w_s_delta = (1.0-b2c2)/denom;
    w_a_delta = 2.0*((1.0-otb2)/denom)*cos_theta_cs_mc;
    w_a_opp = 2.0 * ((1.0+otb2+one_m_b2+alpha*(1.0-otb2))/denom) * (-1.0*cos_theta_cs_mc); 
    w_s_xi_opp = one_m_b2/denom;
    w_a_xi_opp = 2.0*(one_m_b2/denom)*(-1.0*cos_theta_cs_mc);
    w_s_delta_opp = (1.0-b2c2)/denom;
    w_a_delta_opp = 2.0*((1.0-otb2)/denom)*(-1.0*cos_theta_cs_mc);
    # return weights
    tmp_w_list = [ w_a,w_a_opp,w_s_xi,w_s_xi_opp,w_a_xi,w_a_xi_opp ]
    tmp_w_list += [ w_s_delta,w_s_delta_opp,w_a_delta,w_a_delta_opp ]
    tmp_w_list += [beta_mc]
    return tmp_w_list

