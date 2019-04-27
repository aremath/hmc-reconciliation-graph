class Histogram:
    def __init__(self, init):
        '''
        Initialize the Histoogram to either be {} or {0:1}. {0:1} corresponds to (1a) in the paper,
        and {} corresponds to (1b) in the paper.
        '''
        assert(init == None or init == 0 or type(init) == dict)
        if init == 0:
            self.histogram_dict = {0:1}
        elif type(init) == dict:
            self.histogram_dict = init
        else:
            self.histogram_dict = {}
    
    def shift(self, value):
        '''
        Shift the histogram by some value.
        This corresponds to +1 and +2 constants in the paper.
        {0:1, 3:2} << 2   => {2:1, 5:2}
        '''
        new_dict = {}
        old_dict = self.histogram_dict
        for old_key in old_dict :
            new_dict[old_key + value] = old_dict[old_key]
        return Histogram(new_dict)

    def combine(self, other):
        '''
        Add the two histogram together, the way you would expect it to work.
        This corresponds to the max operation in the paper.
        {0:1, 3:2} + {0:3, 5:1}  => {0:4, 3:2, 5:1}
        '''
        new_hist = Histogram.sum([self, other])
        return new_hist
    
    @staticmethod
    def sum(hist_list):
        new_dict = {}
        for hist in hist_list :
            for key in hist.histogram_dict :
                if key not in new_dict :
                    new_dict[key] = 0
                new_dict[key] += hist.histogram_dict[key]
        return Histogram(new_dict)
    
    def product_combine(self, other):
        '''
        {0:1} * {3:1}  => {3:1}
        {0:1, 2:1} * {3:1}  => {3:1. 5:1}
        {0:1, 2:1} * {3:1, 5:2}  => {3:1, 5:3, 7:2}
        #TODO: Update the documentation
               We multiply the value by 2 if it's not
               doing anything with a 0 key.
        '''
        new_dict = {}

        old_dict_A = self.histogram_dict
        old_keys_A = old_dict_A.keys()
        old_dict_B = other.histogram_dict
        old_keys_B = old_dict_B.keys()

        for old_key_A in old_keys_A :
            for old_key_B in old_keys_B :
                new_key = old_key_A + old_key_B
                if new_key not in new_dict:
                    new_dict[new_key] = 0
                if old_key_A == 0 or old_key_B == 0:
                    new_dict[new_key] += old_dict_A[old_key_A] * old_dict_B[old_key_B]
                else:
                    new_dict[new_key] += old_dict_A[old_key_A] * old_dict_B[old_key_B] * 2

        return Histogram(new_dict)
    
    def __eq__(self, other):
        return self.histogram_dict == other.histogram_dict
    
    def __ne__(self, other):
        return not self == other
    
    def __lshift__(self, value):
        return self.shift(value)
    
    def __add__(self, other):
        return self.combine(other)
    
    def __mul__(self, other):
        return self.product_combine(other)

    def __repr__(self):
        return str(self.histogram_dict)

 # Tests

if __name__ == '__main__':

    def testShift1():
        hist = Histogram(None)
        new_hist = hist << 1

        assert(new_hist.histogram_dict == {})
    
    def testShift2():
        hist = Histogram(0)
        new_hist = hist << 3

        assert(hist.histogram_dict == {0:1})
        assert(new_hist.histogram_dict == {3:1})
    
    def testShift3():
        hist = Histogram({0:1, 3:5})
        new_hist = hist << 1

        assert(hist.histogram_dict == {0:1, 3:5})
        assert(new_hist.histogram_dict == {1:1, 4:5})
    
    def testCombine1():
        histA = Histogram(0)
        histB = Histogram(None)
        new_hist = histA + histB

        assert(histA.histogram_dict == {0:1})
        assert(histB.histogram_dict == {})
        assert(new_hist.histogram_dict == {0:1})
    
    def testCombine2():
        histA = Histogram({0:1, 2:1, 3:2})
        histB = Histogram({2:1, 3:1, 10:1})
        new_hist = histA + histB

        assert(histA.histogram_dict == {0:1, 2:1, 3:2})
        assert(histB.histogram_dict == {2:1, 3:1, 10:1})
        assert(new_hist.histogram_dict == {0:1, 2:2, 3:3, 10:1})

    def testProduct1():
        histA = Histogram(0)
        histB = Histogram(None)
        new_hist = histA * histB

        assert(histA.histogram_dict == {0:1})
        assert(histB.histogram_dict == {})
        assert(new_hist.histogram_dict == {})
    
    def testProduct2():
        histA = Histogram({0:1, 2:1, 3:2})
        histB = Histogram({2:1, 3:1, 5:1})
        new_hist = histA * histB

        assert(histA.histogram_dict == {0:1, 2:1, 3:2})
        assert(histB.histogram_dict == {2:1, 3:1, 5:1})
        assert(new_hist.histogram_dict == {2:1, 3:1, 4:1, 5:4, 6:2, 7:1, 8:2})

    testShift1()
    testShift2()
    testShift3()
    testCombine1()
    testCombine2()
    testProduct1()
    testProduct2()
    
