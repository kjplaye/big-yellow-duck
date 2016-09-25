# Need to fix buffer copies in filter mults....numpy.ndarray should talk to ctypes

import numpy
import math
import tempfile
import os
import random
import gnuplot
import exceptions
import ossaudiodev
from ctypes import *

PITCH_NOISE_RATIO = 1.0
PULSE_SHAPE = 1.3

numpy_ptr = numpy.ctypeslib.ndpointer(dtype=numpy.float64,ndim=1,flags='C_CONTIGUOUS')

speech_so = cdll.LoadLibrary('/home/kevin/PYTHON_CODE/_speech.so')

speech_so.zero_filt_mul.restype = None
speech_so.zero_filt_mul.argtypes = [POINTER(c_double),POINTER(c_double),numpy_ptr,c_int,c_int]

speech_so.pole_filt_mul.restype = None
speech_so.pole_filt_mul.argtypes = [POINTER(c_double),POINTER(c_double),numpy_ptr,c_int,c_int]

speech_so.durbin.restype = None
speech_so.durbin.argtypes = [POINTER(c_double),POINTER(c_double)]

speech_so.lpc2cepst.restype = None
speech_so.lpc2cepst.argtypes = [POINTER(c_double),POINTER(c_double)]

speech_so.synth.restype = None
speech_so.synth.argtypes = [POINTER(c_double), c_int, POINTER(c_double), c_double, c_double, c_double, c_double]


def formant_freq(z):
    return 8000 * math.atan2(z.imag,z.real) / (2*math.pi)

def formant_bandwidth(z):
    bw_thresh = 3.0 - 2.0 * 2**0.5
    w = abs(z)
    if w>bw_thresh:
        return 8000 * math.acos(-(w**2 - 4*w + 1)/(2*w)) / math.pi
    else:
        return 8000

def complex_cmp(x,y):
    return cmp(abs(y),abs(x))


class pcm_array(numpy.ndarray):
    """ main class for working with pcm data, an array of floats.
    """
    sample_rate = 8000
    encoding = 'S16_LE'
    def __or__(self,other):
        """ Concatenate
        """
        return pcm(numpy.concatenate((self,other)))
    def autoscale(self):
        """ puts sample between -1 and 1
        """
        m = max(abs(max(self)),abs(min(self)))
        if m:
            return self/m
        else:
            return self
    def strip(self,window = 160,thresh = 0.25):
        """ cuts out energy from begining and end
            we look for a window of samples with energy > thresh * max_energy
        """
        E = self.autoscale()
        E = E * E
        T = numpy.zeros(len(E)-window)
        T[0] = sum(E[0:window])
        for i in xrange(1,len(E)-window):
            T[i] = T[i-1] - E[i-1] + E[i+window-1]
        P = [i for i in xrange(0,len(T)) if T[i] >= thresh]    
        if P:
            return pcm(self[max(0,P[0]-window):min(len(self),P[-1]+2*window)])
    def split(self,window = 160,thresh = 0.25):
        """ 
        """
        E = self.autoscale()
        E = E * E
        T = numpy.zeros(len(E)-window)
        T[0] = sum(E[0:window])
        for i in xrange(1,len(E)-window):
            T[i] = T[i-1] - E[i-1] + E[i+window-1]
        begin_i = 0
        ans = []
        for i in xrange(1,len(T)):
            if T[i] < thresh and T[i-1] >= thresh:
                ans.append(pcm(self[begin_i:i+window]))
                begin_i = None
            elif T[i] >= thresh and T[i-1] < thresh:
                begin_i = i
        return ans

    def reverse(self):
        """ Unlike list.reverse we do not reverse in place.
        """
        x = list(self)
        x.reverse()
        return pcm(x)
    def write_16t(self,filename):
        """ writes a signed 16-bit little endian sample file
        """
        m = max(abs(max(self)),abs(min(self)))
        if m > 1.0:
            data = self.autoscale()
        else:
            data = self
        x = [0 for i in xrange(0,len(data)*2)]
        for i in xrange(0,len(data)):
            val = int(data[i] * 32767)
            if val < 0:
                val+=0x10000
            x[2*i]   = val & 0xff
            x[2*i+1] = (val >> 8) & 0xff
        if isinstance(filename, str):
            f = open(filename,'w')
            f.write(''.join(map(chr,x)))
            f.close()        
        else:
            filename.write(''.join(map(chr,x)))
    def play(self):
        """ play using ossaudiodev
        """
        dsp = ossaudiodev.open('w')
        dsp.write(''.join([chr(int(round(e*126 + 128))) for e in self.autoscale()]))
        dsp.close()
    def play2(self):
        """ play using aplay
        """
        f = tempfile.NamedTemporaryFile()
        self.write_16t(f)
        cmd_line = 'aplay -t raw -f S16_LE -r %d %s 2> /dev/null' % (self.sample_rate,f.name);
        os.system(cmd_line);
        f.close()
    def show(self):
        """ show using wview
        """
        f = tempfile.NamedTemporaryFile()
        self.write_16t(f)
        cmd_line = '/home/kevin/bin/wview %s' % f.name;
        os.system(cmd_line);
        f.close()
    def show2(self):
        """ show using gnuplot
        """
        gnuplot.plot_list(self)
    def show3(self,window = 160):
        """ show gain contour
        """
        E = [e*e for e in self]
        T = [0 for i in xrange(0,len(E)-window)]
        T[0] = sum(E[0:window])
        for i in xrange(1,len(E)-window):
            T[i] = T[i-1] - E[i-1] + E[i+window-1]
        gnuplot.plot_list(T)

def rec(sample_len = 8000,time = None):
    if time:
        to_read = time * 8000
    else:
        to_read = sample_len
    dsp = ossaudiodev.open('r')
    x =dsp.read(to_read)
    dsp.close()
    return pcm([ord(e) - 128 for e in x]).autoscale()        

def pcm(filename=None, format = 'S16_LE'):
    """ pcm(x) loads file x, when x is a string, otherwise uses numpy.array.__init__
        Supported formats:
            S16_LE - signed   16-bit Little-Endian
            S8     - signed    8-bit
            U8     - unsigned  8-bit
             8     - unknown   8-bit
        for unknowns we pick the version which gives the lowest energy after the filter (1-z)
    """
    if isinstance(filename,str):
        f = open(filename,'r')
        data = f.read()
        f.close()
        if format == 'S16_LE':
            x = pcm_array(len(data)/2)
            for i in xrange(0,len(data),2):
                val = ord(data[i]) + ord(data[i+1]) * 256
                if val >= 0x8000:
                    x[i/2] = (val - 0x10000)/32768.0
                else:
                    x[i/2] = val/32768.0
            return x
        elif format == 'S8':
            x = pcm_array(len(data))
            for i in xrange(0,len(data)):
                val = ord(data[i])
                if val >= 0x80:
                    x[i] = (val - 256)/128.0
                else:
                    x[i] = val/128.0
            return x
        elif format == 'U8':
            x = pcm_array(len(data))
            for i in xrange(0,len(data)):
                val = ord(data[i])
                x[i] = (val - 128)/128.0
            return x
        elif format == '8':
            aS = pcm(filename,format = 'S8')
            aS_del = pcm([aS[i+1] - aS[i] for i in xrange(len(aS)-1)])
            aS_del_energy = sum(aS_del * aS_del)

            aU = pcm(filename,format = 'U8')
            aU_del = pcm([aU[i+1] - aU[i] for i in xrange(len(aU)-1)])
            aU_del_energy = sum(aU_del * aU_del)

            return aS if aS_del_energy < aU_del_energy else aU
        else:
            raise Exception('Unsupported format')
            
    else:
        x = pcm_array(len(filename))
        x[:] = filename[:]
        return x


class zero_filt(list):
    def __init__(self,the_coefficients = [1.0]):
        list.__init__(self,the_coefficients)        
        self.memory = [0.0 for i in xrange(0,len(self))]
    def __mul__(self,buffer):
        temp_mem    = (c_double * len(self.memory))()
        temp_coef   = (c_double * len(self))()

        temp_mem[:] = self.memory[:]
        temp_coef[:] = self[:]

        speech_so.zero_filt_mul(temp_mem,temp_coef,buffer,len(self),len(buffer))
        
        self.memory[:] = temp_mem[:]
        self[:] = temp_coef[:]        
        return buffer
    def reset(self):
        self.memory = [0.0 for i in xrange(0,len(self))]

class pole_filt(list):
    def __init__(self,the_coefficients = [1.0]):
        if the_coefficients[0] != 1.0:
            raise Exception('Make sure the first coefficient is 1.0 please')
        list.__init__(self,the_coefficients)
        self.memory = [0.0 for i in xrange(0,len(self))]
    def __mul__(self,buffer):

        temp_mem    = (c_double * len(self.memory))()
        temp_coef   = (c_double * len(self))()

        temp_mem[:] = self.memory[:]
        temp_coef[:] = self[:]

        speech_so.pole_filt_mul(temp_mem,temp_coef,buffer,len(self),len(buffer))

        self.memory[:] = temp_mem[:]
        self[:] = temp_coef[:]                    
        return buffer
    def reset(self):
        self.memory = [0.0 for i in xrange(0,len(self))]


class filt:
    def __init__(self,zero_filter = [1.0], pole_filter = [1.0]):
        self.zero_filter = zero_filt(zero_filter)
        self.pole_filter = pole_filt(pole_filter)
    def __mul__(self,buffer):
        if self.zero_filter != [1.0]:
            self.zero_filter * buffer
        if self.pole_filter != [1.0]:
            self.pole_filter * buffer
        return buffer
    def inverse(self):
        z0 = self.zero_filter[0]
        return filt([e/z0 for e in self.pole_filter],[e/z0 for e in self.zero_filter])
    def show(self):
        """ display the log spectral magnitude for the filter
        """
        dots = 4000
        z = (math.exp(1)) ** (1J * math.pi / dots)
        A = [math.log(abs(sum([z**((i+1)*j) * self.zero_filter[j] for j in xrange(0,len(self.zero_filter))]))) for i in xrange(0,dots)]
        B = [math.log(abs(sum([z**((i+1)*j) * self.pole_filter[j] for j in xrange(0,len(self.pole_filter))]))) for i in xrange(0,dots)]
        gnuplot.plot_list([A[i] - B[i] for i in xrange(0,dots)])
    def reset(self):
        self.zero_filter.reset()
        self.pole_filter.reset()
    def chain_memory(self,x):
        """ chain the memory from a previous filter
        """
        if isinstance(x,filt):
            self.zero_filter.memory = x.zero_filter.memory
            self.pole_filter.memory = x.pole_filter.memory
        else:
            raise Exception('Can only chain filters into filters')

# 2nd order Butterworth 100 Hz HPF:
hpf = filt([0.946, -1.892, 0.946],[1.0, -1.889033, 0.8948743])

class lp_model:
    '''  Linear predictive model for speech

         Applying 1-\sum_{i=1}^order a_i z^i minimizes error
         a_i are the linear predictive coefficients

         Should do formats:
            pcm: Assuming we are seeing a pre-emphisized,
                 filtered cut of speech...we hit it with a Hamming window here
            lpc: Linear prdictive coeffiecients
            rc: Refletion coeffiecients
            ac: Auto correlations
            lsf: Line spectral sequences
            lar: Log area ratios
            roots: Complex roots


         to convert between formats, use
            lp_model(format = YOUR_FORMAT) and
            conversion methods

         some methods
         auto-convert...for instance show creates roots
    '''
    order = 10
    hamming_window = dict() #All instances *should* update this
    def ac2lpc(self):
        #Do Durbin Recursion        
        self.lpc = [0 for i in xrange(0,self.order)]

        temp_ac = (c_double * len(self.ac))()
        temp_lpc = (c_double * self.order)()
        temp_ac[:] = self.ac[:]
        
        speech_so.durbin(temp_ac,temp_lpc)

        self.lpc[:] = temp_lpc[:]

    def lpc2cepst(self):
        self.cepst = [0 for i in xrange(0,self.order)]

        temp_cepst = (c_double * self.order)()
        temp_lpc = (c_double * self.order)()
        temp_lpc[:] = self.lpc[:]

        speech_so.lpc2cepst(temp_lpc,temp_cepst)
        
        self.cepst[:] = temp_cepst[:]
        
    def lpc2roots(self, max_try = 100, max_iter = 100, root_eps = 0.001):
        # Newton's method: 
        #  x ---> x - f(x)/f'(x)
        roots_found = 0;
        self.root = [0 for i in xrange(0,self.order)]
        for i in xrange(0,max_try):
            if roots_found >= self.order:
                break            
            x = ((4.0*random.random()) - 2) + ((4.0*random.random()) - 2)*1J
            for j in xrange(0,max_iter):
                last_x = x

                pow_x = 1
                f = 1
                f_prime = 0
                for k in xrange(0,self.order):
                    f_prime += pow_x * -self.lpc[k] * (k+1);
                    pow_x *= x;
                    f += pow_x * -self.lpc[k];

                x = x - f / f_prime;

                dist = abs(x - last_x)

                if (dist > 1/root_eps):
                    break
                if (dist < root_eps):
                    if (x.imag < 0):
                        x = x.conjugate()
                    no_match_flag = True
                    for k in xrange(0,roots_found):
                        if (abs(x - self.root[k]) < root_eps):
                            no_match_flag = False
                            break
                    if (no_match_flag):
                        self.root[roots_found] = x
                        roots_found += 1
                    break
        self.root = [1.0/e.conjugate() for e in self.root[0:roots_found]]
        self.root.sort(complex_cmp)

    def __repr__(self):
        if hasattr(self,'root'):
            return 'lp_model with formants: %s HZ, bandwidths: %s HZ' % (', '.join(['%4.1f' % formant_freq(e) for e in self.root]), ', '.join(['%4.1f' % formant_bandwidth(e) for e in self.root]))
        else:
            return 'lp_model with LPC = %s' % ', '.join(['%3.3f' % e for e in self.lpc])
        
    def __init__(self, data=None, order=10, format='pcm'):
        self.order = order
        if isinstance(data,pcm_array):
            #Throw down a Hamming Window
            if len(data) not in lp_model.hamming_window:
                lp_model.hamming_window[len(data)] = numpy.hamming(len(data))
            new_data = data * lp_model.hamming_window[len(data)]

            #Find auto correlations
            self.ac = numpy.correlate(numpy.concatenate((new_data,numpy.zeros(self.order))),new_data)
            self.ac2lpc()
            self.lpc2cepst()

    def vocal_tract_filter(self):
        """ returns the model vocal tract filter
        """
        return filt([1.0],[1.0] + [-e for e in self.lpc])

    def show(self):
        """ display the roots of the vocal tract filter
                  the red circle is complex modulus 1
                  the green circle represents a band width of about 200 Hz
        """
        if (not(hasattr(self,'root'))):
            self.lpc2roots()
        R = [[r.real,r.imag] for r in self.root + [e.conjugate() for e in self.root]]
        bw_rad = 0.9245
        dots = 1000
        C1 = [[bw_rad * math.cos(2*3.141*i/dots), bw_rad * math.sin(2*3.141*i/dots)] for i in xrange(0,dots)]
        C2 = [[math.cos(2*3.141*i/dots), math.sin(2*3.141*i/dots)] for i in xrange(0,dots)]
        gnuplot.plot_lists([C2,C1,R],the_options=['title \"\"' for i in xrange(0,3)])

    def show2(self):
        """ display the log spectral magnitude of the vocal tract filter
        """
        self.vocal_tract_filter().show()    

pulse_160 = pcm([math.exp(-PULSE_SHAPE*x)*x for x in range(160)])
pulse_energy = math.sqrt(sum(pulse_160 * pulse_160))

# pulse = pulse with energy a
#     N = noise with energy b
# If X = n pulses + N then
# energy(X) = a*n + b*T
#     ac(X) = a*n
# a = ac / n
# b = (en - ac) / T

class excitation_model:
    def __init__(self, data=None):
        if isinstance(data,pcm_array):
            auto_c = numpy.correlate(numpy.concatenate((data,numpy.zeros(len(data)))),data)

            maxi = 20
            the_max = 0;
            ans = 0
            for i in xrange(20,len(auto_c)):
                temp = auto_c[i]
                if temp > the_max:
                    the_max = temp
                    maxi = i
            for i in xrange(len(data)):
                temp = data[i]
                ans = ans + temp * temp

            self.gain = PITCH_NOISE_RATIO * math.sqrt(abs((ans - auto_c[maxi])/len(data)))
            self.pitch = 8000.0 / maxi
            self.gain_on_pitch = math.sqrt(auto_c[maxi]*maxi/len(data))
            
    def __repr__(self):
        return 'excitation_model with gain = %f, pitch = %f, gain_on_pitch %f' % (self.gain,self.pitch,self.gain_on_pitch)

    def synthesize(self,size = 40,time_since_last_pulse = 0):
        ans = (c_double * size)()
        tslp = c_double(time_since_last_pulse)

        speech_so.synth(ans,size,pointer(tslp),self.gain,8000.0 / self.pitch,self.gain_on_pitch,pulse_energy)
        ans = pcm(ans[:])
        return ans, tslp.value

class vocoder(list):
    """ Basic analysis and synthesis of speech
        modify frame_step after encodeing to change the speed post synthesis        
    """
    def __init__(self, data=None,frame_length = 160,frame_step = 40,residual = None):
        self.frame_length = frame_length
        self.frame_step = frame_step
        if isinstance(data,pcm_array):
            new_data = pcm([e for e in data])
            hpf.reset()
            hpf * new_data

            last_filter = None
            if isinstance(residual,list):
                residual[:] = pcm([0.0 for e in data])[:]
            else:
                residual = pcm([0.0 for e in data])
            list.__init__(self,[])
            for i in xrange(self.frame_length/2,len(new_data) - self.frame_length/2,self.frame_step):
                a = i-self.frame_step/2
                b = i+self.frame_step/2
                a2 = i-self.frame_length/2
                b2 = i+self.frame_length/2
                lp_m = lp_model(new_data[a2:b2])
                f = lp_m.vocal_tract_filter().inverse()
                if last_filter:
                    f.chain_memory(last_filter)

                residual[a:b] = new_data[a:b]
                residual[a:b] = f * residual[a:b]
                self.append([lp_m,None])
            count = 0
            for i in xrange(self.frame_length/2,len(new_data) - self.frame_length/2,self.frame_step):
                a = i-self.frame_step/2
                b = i+self.frame_step/2
                a2 = i-self.frame_length/2
                b2 = i+self.frame_length/2
#                ex_m = excitation_model(residual[a2:b2])
                ex_m = excitation_model(new_data[a2:b2])
                self[count][1] = ex_m
                count += 1
        else:
            list.__init__(self,data)
            
    def synthesize(self):
        synth = pcm([])
        last_filter = None
        tslp = 0
        for frame in self:
            ex,tslp = frame[1].synthesize(size = self.frame_step,time_since_last_pulse = tslp)
            f  = frame[0].vocal_tract_filter()
            if last_filter:
                f.chain_memory(last_filter)
            synth = synth | (f * ex)
            last_filter = f
        return synth
        

# TODO:
#       class hmm
    
