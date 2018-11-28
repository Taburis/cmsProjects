

from ROOT import TCanvas, TPad, TStyle

class anPad (object):
    def __init__(self):
        self.width = 400
        self.height = 400
        self.data = []

    def set_size(self, w, h):
        self.height = h
        self.width = w

    def add(self, h, label = ''):
        self.data.append((h,label)) 

    def est_y_range(self):
        if len(self.data) == 0 : return
        for dat in self.data :
            if hasatrr(self, ymax):
                self.ymax = max(dat[0].GetMaximum(),self.ymax)
            if hasatrr(self, ymin):
                self.ymin = min(dat[0].GetMinimum(),self.ymin)
        scale = 10
        dy = (self.ymax - self.ymin)/scale
        self.ymax = self.ymax + dy
        self.ymin = self.ymin - dy

    def set_y_range(self):
        if hasatrr(self, ymax):
            for dat in self.data:
                dat[0].SetAxisRange(self.ymin, self.ymax, "y")
        else : 
            print("please call: est_y_range() to estimate the range ahead")
            return

    def set_style(self, sty):
        self.style = sty

    def draw(self):
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
        for i in range(n*m):
            self.padlist.append(anPad())
            self.padlist[i].set_size(self.width, self.height)

    def flatten(self, i, j):
        return (i-1)*self.ncol+(j-1)

    def set_all_style(self, sty):
        for pad in self.padlist:
            pad.set_style(sty)

    def add(self, i, j, hist, label = ''):
        self.padlist[self.flatten(i,j)].add(hist, label)

    def get_pad(self, i, j):
        return self.padlist[self.flatten(i,j)]
    
    def draw(self):
        self.cd()
        self.Divide(self.ncol, self.nrow)
        for i in range(self.nrow):
            for j in range(self.ncol):
                print(i,j)
                print(self.flatten(i+1,j+1))
                self.cd(self.flatten(i+1,j+1)+1)
                self.get_pad(i+1,j+1).draw()

    def add_tensor_hist(self, matrix):
        for i in range(self.nrow):
            for j in range(self.ncol):
                self.add(i+1,j+1, matrix[i,j])
       
    @classmethod
    def auto_fit(self, tensors, h = 325, w = 350):
        print(tensors[0].shape)
        if tensors[0].nform > 2:
            print('can only support 2D hists array')
            return 
        else : 
            canvs = anCanvas('canvas_'+tensors[0].name, '', tensors[0].shape[0], tensors[0].shape[1], w, h) 
            for tsor in tensors:
                canvs.add_tensor_hist(tsor)
            return canvs

