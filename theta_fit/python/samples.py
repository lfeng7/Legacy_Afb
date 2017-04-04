# a simple class for finding all information about MC samples
import sys

class sample_info(object):
    """
    A tiny class serves really like a structured data, with sample key, xsec, weight etc
    """
    def __init__(self, raw_info):
        super(sample_info, self).__init__()
        self.raw_info = raw_info
        self.key = raw_info[0]
        self.type = raw_info[1]
        self.ngen = raw_info[2]
        self.xsec = raw_info[3]
        self.weight = raw_info[4]

    def print_out(self):
        print 'key = %s, type = %s, ngen = %i, xsec = %.1f, weight = %.3f'%(self.key,self.type,self.ngen,self.xsec,self.weight)


class samples(object):
    """docstring for samples"""
    def __init__(self, txtfile,verbose=False):
        super(samples, self).__init__()
        self.txtfile = txtfile
        self.info = {}
        self.Lumi = 19700
        self.import_info()
        self.verbose = verbose

    def import_info(self):
        """
        Read the txt file with sample name, cross section and nevts generated
        store information in following a hashtable
        key = partial_file_name, value = [sample_type, nevts_gen, xsec]
        """
        # Get input MC and data files according to txt file
        txt_MC = open(self.txtfile)
        for line in txt_MC:
            items = line.split()
            sample_key = items[1]
            sample_type = items[2]
            nevts_gen = int(items[3])
            xsec = float(items[4])
            weight = self.Lumi*xsec/nevts_gen
            if '#' in items[0] or len(items)<3 : continue
            # TT_CT10_qq               qq      21560109        245.9   
            self.info[sample_key] = sample_info([sample_key,sample_type,nevts_gen,xsec,weight])
        txt_MC.close()
        print '(info) Done samples.import_info'

    def get_sample_info(self,filepath):
        """
        input: a string of sample file path
        output: a sample_info object associate with that sample
        """
        key = [ikey for ikey in self.info if ikey in filepath ] # this is hardcoded due to the existence of T_t anf T_tW
        key.sort()
#        print '(Debug) filepath:',filepath
#        print 'key,',key

        if len(key)==0:
            if self.verbose:
                print 'Cannot find information for file %s in txt file %s! Will stop.'%(filepath,self.txtfile)
            return None
        else:
            key = key[-1]
            info_obj = self.info[key]
            info_obj.print_out()
            return info_obj
        

        
