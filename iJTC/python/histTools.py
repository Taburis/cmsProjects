
from ROOT import *
import math

class tensor_hist(object):
    def __init__(self, name_, shape_):
        self.shape = shape_
        self.name = name_
        self.hist = []
        self.dim = 1
        self.nform = len(shape_)
        for i in shape_: 
            self.dim=self.dim*i

    def do_load_hist(self, name, dim, N, f):
        indx = list(dim)
        for i in range(len(dim)):
            indx[i] = 0

        ptr = len(indx)-1
        for i in range(N):
            key = name
            #make hist name
            for j in indx:
                key = key+'_'+str(j)
            indx[ptr]=indx[ptr]+1
            bits = 0
            #make increment
            for k in range(len(dim)-1, -1, -1):
                num = indx[k]+bits
                bits = math.floor(num/dim[k]) 
                indx[k] = int((num)%dim[k])
            print('loading hist: '+key)
            self.hist.append(f.Get(key))

    def load(self, f):
        return self.do_load_hist(self.name, self.shape, self.dim, f)

    def flatten(self, indx):
        n = 0
        if len(indx) != self.nform:
            print('index dimension not match!')
            return
        return self.get_flatten_index(indx)

    def get_flatten_index(self, index):
        indx = 0
        for i in range(len(index)):
            n =1
            for k in range(i+1, len(index)):
                n= n*self.shape[k]
            indx = index[i]*n+indx
        return indx

    def __getitem__(self, index):
        return self.hist[self.flatten(index)]

    def selfLoop(self, func):
        for h in self.hist:
            func(h)

    def reset_name(self,name_):
        for h in self.hist:
            h.SetName(h.GetName().replace(self.name, name_))
        self.name = name_

    def __add__(self, other):
        res = self.clone(self.name+"_added")
        res.add(other)
        return res
    
    def clone(self, newname):
        h = tensor_hist(newname, self.shape)
        h.hist = []
        for hh in self.hist:
            h.hist.append(hh.Clone(hh.GetName().replace(self.name, newname)))
        return h

    def add(self, other):
        if self.shape != other.shape:
            print('Error: try to add two tensor hist with difference shapes')
            return
        print('adding '+self.name+' + '+other.name)
        for i in range(self.dim):
            self.hist[i].Add(other.hist[i])
  
    def SetLineColor(self, color):
        for i in self.hist: i.SetLineColor(color)

    def SetMarkerSize(self, size):
        for i in self.hist: i.SetMarkerSizeColor(size)

    def SetMarkerStyle(self, num):
        for i in self.hist: i.SetMarkerStyle(num)

