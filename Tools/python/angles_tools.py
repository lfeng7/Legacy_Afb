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

