""" A dictionary class that takes integer keys and numerical values
    and returns zero when a nonexistent key is accessed   """

from numbers import Number

class zdict(dict):
    def __init__(self, *args):
        dict.__init__(self, args)
        
    def __getitem__(self,key):

        if isinstance(key,int):

            if key in dict.keys(self):
                val = dict.__getitem__(self, key)
            else:
                val = 0
        else:
            raise TypeError("Key must be integer type") 

        return val

    def __setitem__(self,key,value):

        if isinstance(key,int):
            if isinstance(value,Number):
                dict.__setitem__(self,key,value)
            else:
                raise TypeError("Value must be a number") 
        else:
            raise TypeError("Key must be integer type") 

          


if __name__ == '__main__':

    z = zdict()
    z[1] = 1
    z[-1] = 1

    print(z)
