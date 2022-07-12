import numpy as np
import os


class AsciiConfig(dict):
    ''' Class to load/write the ascii-based grism configuration files '''



    COMMENTS=('#','%','!',';','$')    # comments in the file

    def __init__(self,conffile):
        
        self.conffile=conffile
        self.confpath=os.path.dirname(self.conffile)

        
        self.orders=[]

        valid=lambda s: (len(s)>0 and not any(map(s.startswith,self.COMMENTS)))

        with open(self.conffile,'r') as fp:
            for line in fp:
                line=line.strip()
                if valid(line):
                    tokens=line.split(' ')
                    key=tokens.pop(0)

                    value=[self.retype(t) for t in tokens if t != '']
                    n=len(value)
                    if n ==0:
                        self.orders.append(key.split('_')[1])
                        value=None
                    elif n==1:
                        value=value[0]
                    else:
                        value=np.array(value)
                    self[key]=value


                    #key=tokens[0].upper()
                    #
                    ## remove empty spaces and retype
                    #value = [self.retype(t) for t in tokens[1:] if t !='']
                    #n=len(value)
                    #
                    #if n==0:
                    #    # this would be a ORDER designation, so record it
                    #    # as a new order
                    #    print(value,tokens)
                    #    self.orders.append(key.split('_')[1])
                    #    value=None
                    #
                    #else:
                    #    # this would be a Key/value pair for a real entry
                    #    if n==1:
                    #        # this is a single element list, so grab the lone value
                    #        value=value[0]
                    #    else:
                    #        # this is for a array, so retype it as a np.array
                    #        value=np.array(value)
                    #self[key]=value


    @staticmethod
    def retype(datum):
        try:
            datum=int(datum)
        except ValueError:
            try:
                datum=float(datum)
            except ValueError:
                pass
        return datum


    def write(self,filename=None):
        if filename is None:
            filename=self.conffile

        with open(filename,'w') as fp:
            for k,v in self.items():
                if v is None:
                    print(k,file=fp)
                else:
                    if isinstance(v,np.ndarray):
                        v=' '.join([str(vv) for vv in v])

                    print(f'{k} {v}',file=fp)


if __name__=='__main__':
    cnf=AsciiConf('/Users/rryan/slitlessutils_config/instruments/WFC3IR/g102.conf')
    cnf.write('t.conf')
