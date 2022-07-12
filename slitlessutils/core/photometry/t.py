class Test:

    def __init__(self,filename):
        print('init')
        self.filename=filename
    
    @classmethod
    def read_csv(cls,filename):
        print('read csv')
        return cls(filename)
    
    @classmethod
    def read(cls,filename):
        print('read')
        obj=cls.read_csv(filename)
        return obj
        
if __name__=='__main__':
    x=Test.read('test')

