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



# Get true cos(theta*) using mc tops p4 and mc q z direction. If no qz direction provided, say, in gg,q1q2 case, use regular cs direction def.
def get_true_angles(mc_t,mc_tbar,qz=0):
    
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
    if bisector_vec.Pz()*qz < -0.1 : # aka bisector and true qz are in opposite direction and qz exists!
        bisector_vec *= -1.0
    #find the CS angle
    cos_theta_cs=top_vec.Dot(bisector_vec);

    return [x_f,ttbar_mass_data,cos_theta_cs]

def Get_beta(m1,m2,mtt):
    return math.sqrt(1-2*(m1*m1+m2*m2)/pow(mtt,2)+pow((m1*m1-m2*m2),2)/pow(mtt,4))


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

