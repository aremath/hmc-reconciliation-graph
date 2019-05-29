import HistogramMain

from pathlib import Path
import numpy as np
import signal
from scipy import stats
import sys
import operator

class TimeoutError(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutError

def sample_hist(hist, n):
    s_v = sum(hist.values())
    # Do not sample if the histogram is smaller than the desired sample size
    if s_v < n:
        # The population for each key
        key_population = [[k]*v for k,v in hist.items()]
        # Flatten that list
        r = reduce(operator.concat, key_population)
        return r
    else:
        k = hist.keys()
        v = hist.values()
        # Convert v to a probability distribution
        p_v = [float(i)/s_v for i in v]
        # would use random.choices in 3.6
        return np.random.choice(k, n, p=p_v)

def find_hists(pathstr, d, t, l, timeout=10):
    p = Path(pathstr)
    all_files = [f for f in p.glob("**/*") if f.is_file()]
    tree_files = [f for f in all_files if f.suffix == ".newick"]
    filenames = []
    histograms = []
    for f in tree_files:
        h = None

        # Time out if it's taking too long to calculate the histogram
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(10)
        try:
            h, _ = HistogramMain.calc_histogram(str(f), d,t,l, False)
        except TimeoutError:
            print("{} timed out".format(f))
        signal.alarm(0)

        if h is not None:
            h_d = h.histogram_dict
            s_v = sum(h_d.values())
            # Minimum sample size for statistical testing to make sense
            if s_v >= 20:
                filenames.append(f)
                histograms.append(h_d)
    return filenames, histograms

def normality(hist_sample):
    # p is 1-probability of rejecting null hypothesis
    # So if p < 0.05 we can say with 95% confidence that the sample is not normal
    _,p = stats.normaltest(hist_sample)
    return p

if __name__ == "__main__":
    names, hists = find_hists(sys.argv[1], 2, 3, 1)
    samples = [sample_hist(h, 10000) for h in hists]
    normalities = [normality(s) for s in samples]
    z = zip(names, hists, samples, normalities)
    z = sorted(z, key=lambda x: x[3])
    print("DATA")
    for e in [(str(i[0]), i[3]) for i in z]:
        print(e)

