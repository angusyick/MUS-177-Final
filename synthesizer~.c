#include "m_pd.h"
#include <stdlib.h>
#include <math.h>
#ifdef NT
#pragma warning( disable : 4244 )
#pragma warning( disable : 4305 )
#endif
#define WAVETABLESIZE 1024
#define PI 3.141592653589793
#define SAMPLERATE 44000

static t_class *synthesizer_class;

typedef struct _synthesizer
{
    t_object x_obj; 	    // obligatory header
    t_float x_f;    	    // frequency
	t_float *wavetable;     // holds the values for the lookup table oscillator
	t_float phase;          // phase
	t_float samplerate;     // rate of samples

    t_float wavefold_gain;  // gain value for wavefold
    t_float offset;         // offset value for wavefold
    t_float toggle_wavefold;// wavefold will be on/off if toggle is 0/1

    t_float cutoff;         // cutoff value for allpass filter
    t_float bandwidth;      // bandwidth value for allpass filter
    t_float toggle_allpass; // allpass filter will be on/off if toggle is 0/1

    t_float limit;          // limit value for clipper
    t_float clip_gain;      // gain value for clipper
    t_float toggle_clip_type;//toggles between hard clip and quadratic soft clip
    t_float toggle_clip;    // clipper will be on/off if toggle is 0/1

    t_float depth;          // depth value for tremolo
    t_float freq;           // frequency value for tremolo sine wave
    t_float toggle_trem;    // tremolo will be on/off if toggle is 0/1

    t_float toggle_wave;    // wave will be sine/square if toggle is 0/1
} t_synthesizer;


// quad interpolation routine
static inline float quad_interpolate(t_synthesizer *x)
{
	int truncphase = (int) (x->phase * WAVETABLESIZE);
	float fr = (x->phase * WAVETABLESIZE) - ((float) truncphase);
	float inm1 = x->wavetable[(truncphase - 1) % WAVETABLESIZE];
	float in   = x->wavetable[(truncphase + 0) % WAVETABLESIZE];
	float inp1 = x->wavetable[(truncphase + 1) % WAVETABLESIZE];
	float inp2 = x->wavetable[(truncphase + 2) % WAVETABLESIZE];

	return in + 0.5 * fr * (inp1 - inm1 +
	 fr * (4.0 * inp1 + 2.0 * inm1 - 5.0 * in - inp2 +
	 fr * (3.0 * (in - inp1) - inm1 + inp2)));
}

// function for wavefolding: If toggle is 0, the previous input is passed through to 
// the output. If toggle is 1, wavefolding is applied to the input. In this case, the 
// input is multiplied by the gain (inlet) which is then added to the offset (inlet)
// to make the ingain. Then, the ingain is plugged into a triangle wave (4 harmonic
// approximation) and outputted as a folded wave.
static inline float wave_folding(float input, float gain, float offset, float toggle){
    
    // handles the case where toggle is not a boolean value
    if(toggle != 0){
        toggle = 1;
    }

    // if toggle is off (= 0), don't add wavefolding
    if(!toggle){
        return input;
    }

    // Note: there aren't any illegal values for gain or offset as far as I can see
    // HOWEVER, high values for gain produce "undesirable" output (although still legal)

    // add gain and offset to input value
    float ingain = (gain * input) + offset;

    // plug input with gain and offset into triangle wave (4 harmonic approximation)
    return cos(0.5 * PI * ingain)
        - 1.0/9.0 * cos(1.5 * PI * ingain)
        + 1.0/25.0 * cos(2.5 * PI * ingain)
        - 1.0/49.0 * cos(3.5 * PI * ingain);
}

// global variables to store the value of allpass filters between sample iterations
double x_0 = 0.0; // input
double x_1 = 0.0; // delayed input
double x_2 = 0.0; // delayed input
double y_0 = 0.0; // allpass output
double y_1 = 0.0; // delayed output
double y_2 = 0.0; // delayed output

// function for allpass filter: If toggle is 0, the previous input is passed through to 
// the output. If toggle is 1, a second order allpass filter is applied to the input. 
// The apf is calculated by storing the values of the previous output and applying it 
// to an equation that puts 2 iterations of this together. Bandwidth and cutoff are
// inlets that shape the allpass filter
static inline float allpass(t_synthesizer *x, float input, float toggle){

    // handles the case where toggle is not a boolean value
    if(toggle != 0){
        toggle = 1;
    }

    // if toggle is off (= 0), don't add apf
    if(!toggle){
        return input;
    }

    // rest of function handles calculations for apf

    float cutoff = x->cutoff/x->samplerate;
    float bandwidth = x->bandwidth/x->samplerate;

    // keep cutoff nonzero and between -0.5 and 0.5 to satisfy calculations later
    if(cutoff == 0){
        cutoff = 0.001;
    }else if(cutoff > 0.5){
        cutoff = 0.5;
    }else if(cutoff < -0.5){
        cutoff = -0.5;
    }

    // keep bandwidth nonzero and between -0.5 and 0.5 to satisfy calculations later
    if(bandwidth == 0){
        bandwidth = 0.001;
    }else if(bandwidth >= 0.5){
        bandwidth = 0.499;
    }else if(bandwidth <= -0.5){
        bandwidth = -0.499;
    }

    // calculate valus for d, tf, and c
    double d = -cos(2.0 * PI * (cutoff));
    double tf = tan(PI * (bandwidth)); 
    double c = (tf - 1.0)/(tf + 1.0); 

    // store values for next function call
    // also, calculate y_0 aka the output
    x_0 = input;   
    y_0 = -c * x_0 + (d - d * c) * x_1 + x_2 - (d - d * c) * y_1 + c * y_2;
    x_2 = x_1;
    x_1 = x_0;
    y_2 = y_1;
    y_1 = y_0;

    return y_0;
}

// function for hard clip: Performs hard clip on given input signal using the 
// given limit value. If signal is not within [-limit, limit], clip the signal 
// so that it is between those values. 
static inline float hard_clip(float input, float limit){
    float output;

    // if it's greater than the limit, hard clip
    if(input > limit){
        output = limit;
    // if it's less than the negative of the limit, hard clip
    }else if(input < -1.0 * limit){
        output = -1.0 * limit;
    // else, just return the input
    }else{
        output = input;
    }
    return output;
}

// function for quadratic_soft_clip: Performs quadratic soft clip on the given 
// input signal using the given limit and gain values. If signal is not within 
// [-1, 1], clip the signal so that it is between those values. If the 
// signal IS within those value, perform softclipping using ax^3 + bx^2 + cx + d 
// where x is the input signal multiplied by the gain and a, b, c, and d are 
// scalar coefficients.
static inline float quadratic_soft_clip(float input, float gain){

    // some parameters
    int samplenumber = 0;
    float a = -0.5f;
    float b = 0.0;
    float c = 1.5f;
    float d = 0.0;
    float output;

    float ingain = gain * input; // get the input
 
    // if it's greater than 1, hard clip
    if(ingain > 1.0)
        output = 1.0;
    // if it's less than -1, hard clip
    else if(ingain < -1.0)
        output = -1.0;
    // else, do the softclipping
    else
        output = a * ingain * ingain * ingain
                 + b * ingain * ingain
                 + c * ingain
                 + d;

    return output;
}

// function for clipper: If toggle is 0, the previous input is passed through to 
// the output. If toggle is 1, a clipper is applied to the input. If toggle_type is 0
// perform hard clip and if toggle_type is 1 perform quadratic soft clip.
static inline float clip(float input, float limit, float gain, float toggle, float toggle_type){

    // handles the case where toggle is not a boolean value
    if(toggle != 0){
        toggle = 1;
    }

    // if toggle is off (= 0), don't add clip
    if(!toggle){
        return input;
    }

    // handles the case where type toggle is not a boolean value
    if(toggle_type != 0){
        toggle_type = 1;
    }

    // if toggle is 0, perform hard clip
    // if toggle is 1, perform quadratic soft clip
    if(!toggle_type){
        return hard_clip(input, limit);
    }
    return quadratic_soft_clip(input, gain);
}

// global variable for the tremolo sine wave
float new_sample = 0;

// function for tremolo: If toggle is 0, the previous input is passed through to 
// the output. If toggle is 1, a tremolo is applied to the input. This works by 
// passing a sine wave with given frequency through the tremolo function 
// 1 - depth * (sin/2 + 0.5) using a given depth value. Depth should be between 0 and 1
static inline float tremolo(float input, float depth, float frequency, float toggle){

    // handles the case where toggle is not a boolean value
    if(toggle != 0){
        toggle = 1;
    }

    // if toggle is off (= 0), don't add tremolo
    if(!toggle){
        return input;
    }

    // depth needs to be in the interval [0,1]
    if(depth < 0){
        depth = 0;
    }else if(depth > 1){
        depth = 1;
    }

    // resets new_sample so that it doesn't reach the float limit
    // also increments for the sine wave
    if(new_sample >= 1000){
        new_sample = 0;
    }else{
        new_sample += 0.001;
    }
    
    // sine wave for the modulator of the tremolo
    float new_sin = sinf(new_sample * frequency / (2 * PI));

    // creates the modulator function for tremolo
    float modulator = 1 - depth * (new_sin/ 2 + 0.5);

    // return the input that has been transformed by the modulator;
    return modulator * input;
}

// function to fill the wave table. If x->toggle_wave is 0, sine wave is put in.
// If it is 1, square wave is put in. The function iterates through the wavetable 
// using its size as an input and fills in each value one by one.
static void fill_wave_table(t_synthesizer *x, float size){
    int i;

    for(i = 0; i < WAVETABLESIZE; i++){

        if(x->toggle_wave == 0){
            // Sine Wave
            *(x->wavetable+i) = sinf(2 * PI * (float)i/size);
        }else if(x->toggle_wave == 1){
            // Square Wave
            if(i <= WAVETABLESIZE / 2) *(x->wavetable+i) = -0.9;
            else *(x->wavetable+i) = 0.9;
        }else if(x->toggle_wave == 2){
            // Triangle Wave
            if(i <= WAVETABLESIZE / 2) *(x->wavetable+i) = 2 * (2 * (float)i/WAVETABLESIZE) - 1;
            else *(x->wavetable+i) = 2 * (2 * (WAVETABLESIZE - (float)i)/WAVETABLESIZE) - 1;;
        }else if(x->toggle_wave == 3){
            // Ramp
            *(x->wavetable+i) = 2 * ((float)i/WAVETABLESIZE) - 1;
        }else if(x->toggle_wave == 4){
            // Sawtooth
            *(x->wavetable+i) = 2 * ((WAVETABLESIZE - (float)i)/WAVETABLESIZE) - 1;
        }
    }
}

/* this is the actual performance routine which acts on the samples.
    It's called with a single pointer "w" which is our location in the
    DSP call list.  We return a new "w" which will point to the next item
    after us.  Meanwhile, w[0] is just a pointer to dsp-perform itself
    (no use to us), w[1] and w[2] are the input and output vector locations,
    and w[3] is the number of points to calculate. */
int toggle = 0;

static t_int *synthesizer_perform(t_int *w)
{
	t_synthesizer *x = (t_synthesizer *)(w[1]);
    t_float *freq = (t_float *)(w[2]);
    t_float *out = (t_float *)(w[3]);
    int n = (int)(w[4]);

	// count from 0
	int blocksize = n;
	int i, sample = 0;
	float phaseincrement;
    float wave_fold_in;
    float allpass_in;
    float clip_in;
    float tremolo_in;

    double in, out1 = 0.0;

    // iterates by sample to create the wave
    while (n--)
    {
        if(toggle != x->toggle_wave){
            toggle = x->toggle_wave;
            fill_wave_table(x, (float)WAVETABLESIZE);
        }

		// calculate the phase increment from the frequency
		// and sample rate - this is the number of cycles per sample
		// freq = cyc/sec, sr = samp/sec, phaseinc = cyc/samp = freq/sr
		phaseincrement = *(freq+sample)/x->samplerate;
		
		// increment the phase and make sure it doesn't go over 1.0
		x->phase += phaseincrement;
		while(x->phase >= 1.0f)
			x->phase -= 1.0f;
		while(x->phase < 0.0f)
			x->phase += 1.0f;
        
        // get value from interpolation, pass that through the wavefolding function which
        // returns another value. Then pass the new value through the allpass filter function
        // and output it.
        wave_fold_in = quad_interpolate(x);
        allpass_in = wave_folding(wave_fold_in, x->wavefold_gain, x->offset, x->toggle_wavefold);
        clip_in = allpass(x, allpass_in, x->toggle_allpass);
        tremolo_in = clip(clip_in, x->limit, x->clip_gain, x->toggle_clip, x->toggle_clip_type);
        *(out+sample) = tremolo(tremolo_in, x->depth, x->freq, x->toggle_trem);

        // increment the sample
		sample++;
    }
    return (w+5);
}

// function to display oscillator in the external
static void synthesizer_dsp(t_synthesizer *x, t_signal **sp)
{
	// we'll initialize samplerate when starting up
	x->samplerate = sp[0]->s_sr;
    dsp_add(synthesizer_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n);
}

static void *synthesizer_new(void)
{
	float size;
    t_synthesizer *x = (t_synthesizer *)pd_new(synthesizer_class);
    outlet_new(&x->x_obj, gensym("signal"));
	// initialize variables
    x->x_f = 0.0f;
	x->phase = 0.0f;
	size = (float)WAVETABLESIZE;

    // create inlets
    // these 3 are for wavefolding
    floatinlet_new(&x->x_obj, &x->wavefold_gain);
    floatinlet_new(&x->x_obj, &x->offset);
    floatinlet_new(&x->x_obj, &x->toggle_wavefold);

    // these three are for filter
    floatinlet_new(&x->x_obj, &x->cutoff);
    floatinlet_new(&x->x_obj, &x->bandwidth);
    floatinlet_new(&x->x_obj, &x->toggle_allpass);

    // these four are for clipping
    floatinlet_new(&x->x_obj, &x->limit);
    floatinlet_new(&x->x_obj, &x->clip_gain);
    floatinlet_new(&x->x_obj, &x->toggle_clip_type);
    floatinlet_new(&x->x_obj, &x->toggle_clip);

    // these three are for tremolo
    floatinlet_new(&x->x_obj, &x->depth);
    floatinlet_new(&x->x_obj, &x->freq);
    floatinlet_new(&x->x_obj, &x->toggle_trem);

    // this one is to change the wave
    floatinlet_new(&x->x_obj, &x->toggle_wave);
	
	// space for WAVETABLESIZE samples
	x->wavetable = (t_float *)malloc(WAVETABLESIZE * sizeof(t_float));
	
    // fills the wave table with values for either sine or square wave
    fill_wave_table(x, size);
	
    return (x);
}

// delete allocated memory
static void synthesizer_free(t_synthesizer *x)
{
	free(x->wavetable);
}

// setup external
void synthesizer_tilde_setup(void)
{
    synthesizer_class = class_new(gensym("synthesizer~"), (t_newmethod)synthesizer_new, (t_method)synthesizer_free,
    	sizeof(t_synthesizer), 0, A_DEFFLOAT, 0);
    CLASS_MAINSIGNALIN(synthesizer_class, t_synthesizer, x_f);
    class_addmethod(synthesizer_class, (t_method)synthesizer_dsp, gensym("dsp"), 0);
}
