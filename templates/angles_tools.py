# Make cos(theta*) and w_a,w_s from the reco t and tbar p4
from ROOT import *
from math import sqrt

# Get cos(theta*) for reco top pair p4
def get_angles(reco_t,reco_tbar):
    
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
    top_data = TVector3(Top_Data.Px(),Top_Data.Py(),Top_Data.Pz())
    proton_data = TVector3(Proton_data.Px(),Proton_data.Py(),Proton_data.Pz())
    proton2_data = TVector3(Proton2_data.Px(),Proton2_data.Py(),Proton2_data.Pz())
    
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

    return [x_f,ttbar_mass_data,cos_theta_cs_data]

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
    delta_mag = proton_data.Mag() - proton2_data.Mag()
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
    quark = q_p4.Clone()
    antiquark = qbar_p4.Clone()
    
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
    
    #Doing the boost
    R_data = R_data.Boost(Bx_data,By_data,Bz_data);
    Top_Data = R_data*Top_Data;
    ATop_Data = R_data*ATop_Data;
    quark = R_data*quark;
    antiquark = R_data*antiquark;
    #Reset the boost
    R_data=S;
    #Define three vectors for P,Pbar,top, quark and antiquark in ttbar c.m frame
    top_data = TVector3(Top_Data.Px(),Top_Data.Py(),Top_Data.Pz())
    true_quark_direction = TVector3(quark.Px(),quark.Py(),quark.Pz())
    true_antiquark_direction = TVector3(antiquark.Px(),antiquark.Py(),antiquark.Pz())
    
    true_antiquark_direction = -1.0*true_antiquark_direction;
    
    #Normalize vectors
    top_data = top_data*(1.0/top_data.Mag());
    true_quark_direction = true_quark_direction*(1.0/true_quark_direction.Mag());
    true_antiquark_direction = true_antiquark_direction*(1.0/true_antiquark_direction.Mag());
    #find the unit bisectors
    bisector_mc = (true_quark_direction+true_antiquark_direction)*(1.0/(true_quark_direction+true_antiquark_direction).Mag());
    #find the CS angle
    cos_theta_cs_true=top_data.Dot(bisector_mc);

    return [x_f,ttbar_mass_data,cos_theta_cs_true]

# Get cos(theta*) for top pair p4 using mc p4 and mc init particles
# Use quark direction as positive direction
def get_true_angles_v2(reco_t,reco_tbar,q_pz,qbar_pz,q_id,qbar_id):
    
    #Lorentz Declartions

    sqrt_s=8000;
    beam_energy=sqrt_s/2;

    S = TLorentzRotation();

    #intialize the rotations to the identity
    R_data = TLorentzRotation() ;
      
    #Vectors of beam1 and beam2 for unit bisector purpose
    # TLorentzVector(px,py,pz,E)
    Top_Data = TLorentzVector(0,0,1,2);
    ATop_Data = TLorentzVector() ;

    # Set up t and tbar p4
    Top_Data = reco_t.Clone()
    ATop_Data = reco_tbar.Clone()

    #assign the quark and antiquark fourvectors from the last two particles in the list of mc_ particles
    quark = (q_pz,q_id)
    antiquark = (qbar_pz,qbar_id)
    
    #Make the 4-vector of the ttbar
    Q_Data = Top_Data + ATop_Data;
    ttbar_mass_data=Q_Data.Mag(); # to return Mtt

    #defining the Px, Py,and Pz, and energies to boost into the ttbar rest frame
    Bx_data = -1*Q_Data.Px()/Q_Data.E();  
    By_data = -1*Q_Data.Py()/Q_Data.E();  
    Bz_data = -1*Q_Data.Pz()/Q_Data.E();
    Qt = sqrt(Q_Data.Px()*Q_Data.Px()+Q_Data.Py()*Q_Data.Py());

    #Feynman x
    x_f = (quark.E()-antiquark.E())/beam_energy
    
    #Doing the boost
    R_data = R_data.Boost(Bx_data,By_data,Bz_data);
    Top_Data = R_data*Top_Data;
    ATop_Data = R_data*ATop_Data;
    #Reset the boost
    R_data=S.Clone();
    #Define three vectors for P,Pbar,top, quark and antiquark in ttbar c.m frame
    top_data = TVector3(Top_Data.Px(),Top_Data.Py(),Top_Data.Pz())
        
    #Normalize vectors
    top_data = top_data*(1.0/top_data.Mag());
    true_quark_direction = true_quark_direction*(1.0/true_quark_direction.Mag());
    #find the CS angle
    cos_theta_cs_true=top_data.Dot(true_quark_direction);

    return [x_f,ttbar_mass_data,cos_theta_cs_true]

# Get reweighting weights for qqbar->ttbar templates
def GetAnglesWeights(Top_MC,ATop_MC,cos_theta_cs_mc,alpha=-0.129):
    M2_1_mc   = Top_MC.Mag2();
    M2_2_mc   = ATop_MC.Mag2();
    ttbar_mass_mc = (Top_MC+ATop_MC).M()
    num_mc    = 1. - 2.*(M2_1_mc+M2_2_mc)/(ttbar_mass_mc*ttbar_mass_mc) + (M2_1_mc-M2_2_mc)*(M2_1_mc-M2_2_mc)/(ttbar_mass_mc*ttbar_mass_mc*ttbar_mass_mc*ttbar_mass_mc);
    denom_1_mc   = (1. + (M2_1_mc-M2_2_mc)/(ttbar_mass_mc*ttbar_mass_mc))*(1. + (M2_1_mc-M2_2_mc)/(ttbar_mass_mc*ttbar_mass_mc));
    denom_2_mc   = (1. + (M2_2_mc-M2_1_mc)/(ttbar_mass_mc*ttbar_mass_mc))*(1. + (M2_2_mc-M2_1_mc)/(ttbar_mass_mc*ttbar_mass_mc));
    beta_mc   = sqrt(sqrt((num_mc*num_mc)/(denom_1_mc*denom_1_mc) * (num_mc*num_mc)/(denom_2_mc*denom_2_mc)));
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
    alpha = -0.129;
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

