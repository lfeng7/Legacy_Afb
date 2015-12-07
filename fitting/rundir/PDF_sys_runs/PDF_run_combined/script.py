import os
from ROOT import *

names = ['Powheg_qq_semilep_TT','Powheg_gg_semilep_TT','Powheg_dilep_TT','Powheg_had_TT','Powheg_qq_semilep_TT_SC','Powheg_gg_semilep_TT_SC','Powheg_dilep_TT_SC','Powheg_had_TT_SC','DY1Jets','DY2Jets','DY3Jets','DY4Jets','W1Jets','W2Jets','W3Jets','W4Jets','T_s-channel','T_t-channel','T_tW-channel','Tbar_s-channel','Tbar_t-channel','Tbar_tW-channel']

weights = [(245.8/25523595.),(245.8/25523595.),(245.8/25523595.),(245.8/25523595.),(245.8/25523595.),(245.8/25523595.),(245.8/25523595.),(245.8/25523595.),(660.6/23994669.),(215.1/2345857.),(65.79/10655325.),(28.59/5843425.),(6662.8/23038253.),(2159.2/33993463.),(640.4/15507852.),(264.0/13326400.),(3.79/259176.),(56.4/3748155.),(11.1/495559.),(1.76/139604.),(30.7/1930185.),(11.1/491463.)]

for i in range(len(names)) :
    name = names[i]
    weight = weights[i]
    file = TFile(name+'_all_histos.root')
    file_sb = TFile(name+'_sb_all_histos.root')
    p_4j = (file.Get(name+'_all_s_p_4j')).Integral()
    p_5j = (file.Get(name+'_all_s_p_5j')).Integral()
    m_4j = (file.Get(name+'_all_s_m_4j')).Integral()
    m_5j = (file.Get(name+'_all_s_m_5j')).Integral()
    p_4j_sb = (file_sb.Get(name+'_sb_all_s_p_4j')).Integral()
    p_5j_sb = (file_sb.Get(name+'_sb_all_s_p_5j')).Integral()
    m_4j_sb = (file_sb.Get(name+'_sb_all_s_m_4j')).Integral()
    m_5j_sb = (file_sb.Get(name+'_sb_all_s_m_5j')).Integral()
    print name+' :'
    print '     L+ :  '+str((25523595./245.8)*weight*(p_4j+p_5j))
    print '     L- :  '+str((25523595./245.8)*weight*(m_4j+m_5j))
    print '     5J :  '+str((25523595./245.8)*weight*(p_5j+m_5j))
    print '     4J :  '+str((25523595./245.8)*weight*(p_4j+m_4j))
    print name+'_sb :'
    print '     L+ :  '+str((25523595./245.8)*weight*(p_4j_sb+p_5j_sb))
    print '     L- :  '+str((25523595./245.8)*weight*(m_4j_sb+m_5j_sb))
    print '     5J :  '+str((25523595./245.8)*weight*(p_5j_sb+m_5j_sb))
    print '     4J :  '+str((25523595./245.8)*weight*(p_4j_sb+m_4j_sb))
