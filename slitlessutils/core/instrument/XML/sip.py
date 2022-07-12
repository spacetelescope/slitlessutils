

class SIP(dict):
    ''' class to hold the SIP data '''

    COMP = ('a','b','ap','bp')
    
    def __init__(self):
        for comp in self.COMP:
            self[comp]={}

    @classmethod
    def from_xml(cls,xml):
        obj=cls()
        for coef in xml:
            obj.append(coef.tag,coef.attrib['i'],coef.attrib['j'],coef.text)
        return obj
            
    def append(self,c,i,j,v):
        if c in self.COMP:
            self[c][(int(i),int(j))]=float(v)
        
    def update_header(self,hdr):
        for comp in self.COMP:
            data=self[comp]
            if data:
                if hdr['CTYPE1'][-3:] != 'SIP':
                    hdr['CTYPE1']+='-SIP'
                if hdr['CTYPE2'][-3:] != 'SIP':
                    hdr['CTYPE2']+='-SIP'
                
                hdr[f'{comp}_order']=(0,'SIP Order')
                order=0
                for (i,j),v in data.items():
                    hdr[f'{comp}_{i}_{j}']=v
                    order=max(max(i,j),order)
                hdr[f'{comp}_order']=order

                

if __name__=='__main__':
    s=SIP()
    s.append('a',2,3,2.)
    s.append('b',1,3,2.)
    hdr={}
    s.update_header(hdr)
    print(hdr)
