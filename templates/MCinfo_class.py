# define a class that read in MC information from a txt file, store each entry in another class, and provide method to 
# find any entry with a filename such as DY4JetsToLL_M

import sys

class MC_info:
    """docstring for ClassName"""

    def __init__(self,txtFilePath):
        # Definition of MC file class
        self.txtFilePath = txtFilePath
        # attributes
        self.chart = {}
        # methods
        self.import_txt()
        # maybe others

    def import_txt(self):
        # read a txt file and parse through
        txt_MC = open(self.txtFilePath)
        for line in txt_MC:
            items = line.split()
            if items[0] == '#' or len(items)<3:
                continue
            # file_name                 key               type      nevts_gen      xsec             title
            # TT_CT10_qq.root         TT_CT10_qq          qq      21560109        245.9        q#bar{q} #rightarrow t#bar{t}
            items = items[1:]
            sample_key = items.pop(0)
            # create a new MC_entry object as a container for current entry
            entry = self.MC_entry(sample_key)
            entry.sample_key = sample_key
            entry.type = items.pop(0)
            entry.nevts_gen = items.pop(0)
            entry.xsec = items.pop(0)
            entry.title = ' '.join(items)
            self.chart[sample_key] = entry
        print 'Importing MC info from %s done!'%self.txtFilePath
    
    def get_entry(self,sample_name):
        # use a key to find the entry stored in import_txt
        entry = [ self.chart[key] for key in self.chart if key in sample_name] # for loop in chart will loop over keys in dic
        if len(entry)==1:
            return entry[0]
        elif len(entry)==0:
            print 'Find no entry for sample %s!'%sample_name
            return None
        else:
            print 'Find %i entry for sample %s!'%(len(entry),sample_name)
            matched_keys = [key for key in self.chart if key in sample_name]
            print 'keys: %s'%' '.join(matched_keys)
            correct_key = max(matched_keys)
            print 'Will return entry of key %s'%correct_key
            return self.chart[correct_key]

    def print_chart(self):
        # print a nicely formatted chart 
        all_sample_keys = self.chart.keys()
        all_sample_keys.sort()
        for ikey in all_sample_keys:
            value = self.chart[ikey]
            print '%s\t\t%s\t%s\t%s\t%s'%(value.sample_key, value.type, value.nevts_gen, value.xsec, value.title)


    # basically a struc to store each entry in MC.txt
    class MC_entry:

        def __init__(self,sample_key):
            self.sample_key = sample_key

if __name__ == '__main__':

    argv = sys.argv[1:]

    # help message
    if len(argv)==0:
        print """
        A class that read and store MC.txt 
        Usage: ./MCinfo_class.py  path_to_MC_txt_file
        """
        sys.exit(1)    

    txt_file = argv.pop(0)
    self = MC_info(txt_file)


        
