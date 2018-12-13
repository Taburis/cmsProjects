

from ROOT import TCanvas, TPad, TStyle

class anPad (object):
    def __init__(self):
        self.width = 400
        self.height = 400
        self.data = []
        self.control = {'yrange':'auto'}
        self.ymax = 0
        self.ymin = 0
        '''
        yrange: auto (automatic setup the yrange) 
        '''

    def set_size(self, w, h):
        self.height = h
        self.width = w

    def add(self, h, label = ''):
        self.data.append((h,label)) 

    def est_y_range(self):
        if len(self.data) == 0 : return
        for dat in self.data :
                self.ymax = max(dat[0].GetMaximum(),self.ymax)
                self.ymin = min(dat[0].GetMinimum(),self.ymin)
        scale = 10
        dy = (self.ymax - self.ymin)/scale
        self.ymax = self.ymax + dy
        self.ymin = self.ymin - dy

    def auto_set_y_range(self):
        for dat in self.data:
            dat[0].SetAxisRange(self.ymin, self.ymax, "y")

    def range_control(self, pset):
        if pset == 'auto' : 
            self.est_y_range()
            self.auto_set_y_range()

    def set_style(self, sty):
        self.style = sty

    def draw(self):
        self.range_control(self.control['yrange'])
        opt = ""
        for h in self.data:
            print(h[0].GetName())
            h[0].Draw(opt)
            opt = "same"

class anCanvas(TCanvas):
    def __init__(self, name, title, n=1, m= 1, h= 325, w=350):
        super(anCanvas, self).__init__(name, title, h*m, w*n)
        self.nrow = n
        self.ncol = m
        self.width = w
        self.height = h
        self.padlist = []
        self.mapper = self.identity_mapper
        for i in range(n*m):
            self.padlist.append(anPad())
            self.padlist[i].set_size(self.width, self.height)

    def identity_mapper(self, i, j):
        return (i, j)

    def remapper(self, i, j):
        x = self.mapper(i,j)
        return self.flatten_mapper(x[0], x[1])

    def set_mapper(self, f):
        self.mapper = f

    def flatten_mapper(self, i, j):
        return (i-1)*self.ncol+(j-1)

    def set_all_style(self, sty):
        for pad in self.padlist:
            pad.set_style(sty)

    def add(self, i, j, hist, label = ''):
        #print("get pad ("+str(i)+", "+str(j)+") : "+str(self.remapper(i,j)))
        self.padlist[self.remapper(i,j)].add(hist, label)

    def get_pad(self, i, j):
        return self.padlist[self.remapper(i,j)]
    
    def draw(self):
        self.cd()
        self.Divide(self.ncol, self.nrow)
        for i in range(self.nrow):
            for j in range(self.ncol):
                self.cd(self.remapper(i+1,j+1)+1)
                self.get_pad(i+1,j+1).draw()

    def add_tensor_hist(self, matrix, f = 0):
        if f == 0 : f = self.identity_mapper
        print(matrix.name)
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                x = f(i+1, j+1)
                self.add(x[0], x[1], matrix[i,j])

    def update_common_labels(self, labels, x, tex):
        for i in range(self.nrow):
            for j in range(self.ncol):
                self.cd(self.remapper(i+1,j+1)+1)
                tex.DrawLatexNDC(x[0], x[1], labels[i][j])

       
    @classmethod
    def auto_fit(self, tensors, shape = 0, fmap = 0, h = 325, w = 350):
        if shape == 0: shape = tensors[0].shape
        print(tensors[0].shape)
        if tensors[0].nform > 2:
            print('can only support 2D hists array')
            return 
        else : 
            canvs = anCanvas('canvas_'+tensors[0].name, '', shape[0], shape[1], h, w) 
            for tsor in tensors:
                canvs.add_tensor_hist(tsor, fmap)
            return canvs

