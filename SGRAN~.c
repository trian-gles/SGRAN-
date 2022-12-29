/**
	@file
	sgran~ - a granulation algorithm designed by Mara Helmuth, rewritten and ported to Max by Kieran McAuliffe

*/

#include "ext.h"
#include "z_dsp.h"
#include "math.h"
#include "ext_buffer.h"
#include "ext_atomic.h"
#include "ext_obex.h"
#include <stdlib.h>
#include <stdbool.h>

#define MAXGRAINS 1000
#define MIDC_OFFSET (261.62556530059868 / 256.0)
#define M_LN2	0.69314718055994529



typedef struct Grain {
	float waveSampInc; 
	float ampSampInc; 
	float wavePhase; 
	float ampPhase; 
	int dur; 
	float panR; 
	float panL; 
	int currTime; 
	bool isplaying;
	} Grain;

typedef struct _sgran {
	t_pxobject w_obj;
	t_buffer_ref *w_buf;
	t_buffer_ref *w_env;
	t_symbol *w_name;
	t_symbol *w_envname;
	
	t_bool running;
	Grain grains[MAXGRAINS]
	
	long w_len;
	long w_envlen;
	
	double freqLow;
	double freqMid;
	double freqHigh;
	double freqTight;

	double grainDurLow;
	double grainDurMid;
	double grainDurHigh;
	double grainDurTight;

	double panLow;
	double panMid;
	double panHigh;
	double panTight;
	
	double grainRateVarLow;
	double grainRateVarMid;
	double grainRateVarHigh;
	double grainRateVarTight;
	
	int newGrainCounter;
	float grainRate;
	
	short w_connected[2];
} t_sgran;




void *sgran_new(t_symbol *s,  long argc, t_atom *argv);
void sgran_free(t_sgran *x);
t_max_err sgran_notify(t_sgran *x, t_symbol *s, t_symbol *msg, void *sender, void *data);
void sgran_assist(t_simpwave *x, void *b, long m, long a, char *s);
void sgran_start(t_sgran *x);
void sgran_stop(t_sgran *x);
void sgran_perform64(t_sgran *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam);
void sgran_dsp64(t_sgran *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags);

void sgran_grainrate(t_sgran *x, double rl, double rm, double rh, double rt);
void sgran_graindur(t_sgran *x, double dl, double dm, double dh, double dt);
void sgran_freq(t_sgran *x, double fl, double fm, double fh, double ft);
void sgran_pan(t_sgran *x, double pl, double pm, double ph, double pt); 

double rrand() 
{
	double min = -1;
	double max = 1;
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

// Taken from RTCMIX code
double octcps(double cps)
{
	return log(cps / MIDC_OFFSET) / M_LN2;
}

double cpsoct(double oct)
{
	return pow(2.0, oct) * MIDC_OFFSET;
}



// From Mara Helmuth
double prob(double low,double mid,double high,double tight)
        // Returns a value within a range close to a preferred value
                    // tightness: 0 max away from mid
                     //               1 even distribution
                      //              2+amount closeness to mid
                      //              no negative allowed
{
	double range, num, sign;

	range = (high-mid) > (mid-low) ? high-mid : mid-low;
	do {
	  	if (rrand() > 0.)
			sign = 1.;
		else  sign = -1.;
	  	num = mid + sign*(pow((rrand()+1.)*.5,tight)*range);
	} while(num < low || num > high);
	return(num);
}

static t_symbol *ps_buffer_modified;
static t_class *s_sgran_class;


void ext_main(void *r)
{
	t_class *c = class_new("sgran~", (method)sgran_new, (method)sgran_free, sizeof(t_sgran), NULL, A_GIMME, 0);

	class_addmethod(c, (method)sgran_dsp64,		"dsp64",	A_CANT, 0);
	class_addmethod(c, (method)sgran_start,		gensym("start"), 0);
	class_addmethod(c, (method)sgran_stop,		gensym("stop"), 0);
	
	class_addmethod(c, (method)sgran_notify,		"notify",	A_CANT, 0);
	
	class_addmethod(c, (method)simpwave_assist,		"assist",	A_CANT, 0);
	
	class_addmethod(c, (method)sgran_grainrate, gensym("grainrate"), 
	A_FLOAT, A_FLOAT, A_FLOAT, A_FLOAT, 0);
	
	class_addmethod(c, (method)sgran_graindur, gensym("graindur"), 
	A_FLOAT, A_FLOAT, A_FLOAT, A_FLOAT, 0);
	
	class_addmethod(c, (method)sgran_freq, gensym("freq"), 
	A_FLOAT, A_FLOAT, A_FLOAT, A_FLOAT, 0);
	
	class_addmethod(c, (method)sgran_pan, gensym("pan"), 
	A_FLOAT, A_FLOAT, A_FLOAT, A_FLOAT, 0);

	class_dspinit(c);
	class_register(CLASS_BOX, c);
	s_sgran_class = c;

	ps_buffer_modified = gensym("buffer_modified");
}
/*
	Inlets:
	0 : On/off
	1 : grain rate
	2 : grain dur
	3 : freq
	4 : pan
	
	1-4 take a list of values for low, mid, high, tight
	
*/

/* Args:
		p0: wavetable
		p1: grainEnv
	*/
void *sgran_new(t_symbol *s,  long argc, t_atom *argv)
{
	t_sgran *x = (t_sgran *)object_alloc(s_sgran_class);
	t_symbol *buf=0;
	t_symbol *env=0;

	dsp_setup((t_pxobject *)x,3);
	buf = atom_getsymarg(0,argc,argv);
	env = atom_getsymarg(1,argc,argv);

	x->w_name = buf;
	x->w_envname = env;
	
	
	//outlets
	outlet_new((t_object *)x, "signal");		// audio outlet l
	outlet_new((t_object *)x, "signal");		// audio outlet r
	
	
	// Setup Grains
	for (size_t i = 0; i < 2000; i++){
        x->grains[i] = (Grain){.waveSampInc=0, 
        	.ampSampInc=0, 
        	.wavePhase=0, 
        	.ampPhase=0, 
        	.dur=0, 
        	.panR=0, 
        	.panL=0, 
        	.currTime=0, 
        	.isplaying=false };
    }

	// create a new buffer reference, initially referencing a buffer with the provided name
	x->w_buf = buffer_ref_new((t_object *)x, x->w_name);
	x->w_env = buffer_ref_new((t_object *)x, x->w_envname);
	
	x->w_len = buffer_getframecount(x->w_buf);
	x->w_envlen = buffer_getframecount(x->w_env);

	return (x);
}




void sgran_free(t_sgran *x)
{
	dsp_free((t_pxobject *)x);

	// must free our buffer reference when we will no longer use it
	object_free(x->w_buf);
	object_free(x->w_env);
}


// A notify method is required for our buffer reference
// This handles notifications when the buffer appears, disappears, or is modified.
t_max_err sgran_notify(t_sgran *x, t_symbol *s, t_symbol *msg, void *sender, void *data)
{
	if (msg == ps_buffer_modified)
		x->w_buffer_modified = true;
	return buffer_ref_notify(x->w_buf, s, msg, sender, data);
}


void sgran_assist(t_sgran *x, void *b, long m, long a, char *s)
{
	if (m == ASSIST_INLET) {	// inlets
		switch (a) {
		case 0:	snprintf_zero(s, 256, "Various messages");	break;
		}
	}
	else if (m == ASSIST_OUTLET){	// outlet
		switch (a) {
		case 0:	snprintf_zero(s, 256, "(signal) right output");	break;
		case 1:	snprintf_zero(s, 256, "(signal) left output");	break;
		}
	}
}


////
// START AND STOP MSGS
////
void sgran_start(t_sgran *x){
	x->running = true;
}

void sgran_stop(t_sgran *x){
	x->running = false;
}


////
// PARAMETER MESSAGES
////

void sgran_grainrate(t_sgran *x, double rl, double rm, double rh, double rt){
	x->grainRateVarLow = rl;
	x->grainRateVarMid = max(rm, rl);
	x->grainRateVarHigh = max(rh, rm);
	x->grainRateVarTight = rt;
}

void sgran_graindur(t_sgran *x, double dl, double dm, double dh, double dt){
	x->grainDurLow = dl;
	x->grainDurMid = max(dm, dl);
	x->grainDurHigh = max(dh, dm);
	x->grainDurTight = dt;
}

void sgran_freq(t_sgran *x, double fl, double fm, double fh, double ft){
	x->freqLow = fl;
	x->freqMid = max(fm, fl);
	x->freqHigh = max(fh, fm);
	x->freqTight = ft;
}

void sgran_pan(t_sgran *x, double pl, double pm, double ph, double pt) {
	x->panLow = pl;
	x->panMid = max(pm, pl);
	x->panHigh = max(ph, pm);
	x->panTight = pt;
}

void new_grain(t_sgran *x, Grain *g){
	
}

void reset_grain_rate(t_sgran *x){
	x->newGrainCounter = (int)round(sys_getsr() * prob(x->grainRateVarLow, x->grainRateVarMid, x->grainDurHigh, x->grainDurTight));
}


// rewrite
void sgran_perform64(t_sgran *x, t_object *dsp64, double **ins, long numins, double **outs, long numouts, long sampleframes, long flags, void *userparam)
{
	t_double		*r_out = outs[0];
	t_double		*l_out = outs[1];
	
	
	int				n = sampleframes;
	float			*b;
	long			len, dex;
	double			v;
	t_buffer_obj	*buffer = buffer_ref_getobject(x->w_buf);
	t_buffer_obj	*env = buffer_ref_getobject(x->w_env);
	t_atom_long		channelcount;

	b = buffer_locksamples(buffer);
	if (!b)
		goto zero;

	channelcount = buffer_getchannelcount(buffer);

	if (x->w_buffer_modified) {
		x->w_buffer_modified = false;
		sgran_limits(x);
		if (!x->w_connected[0])
			min = x->w_start;
		if (!x->w_connected[1])
			max = x->w_end;
	}

	if (min != x->w_start || max != x->w_end) {
		if (min < 0.)
			min = 0.;
		if (max < 0.)
			max = 0.;
		if (min > max)
			x->w_end = x->w_start = min;
		else {
			x->w_start = min;
			x->w_end = max;
		}
		sgran_limits(x);
	}

	b += x->w_begin;
	len = x->w_len;

	if (channelcount == 1) {
		while (n--) {
			v = *in++;
			if (v < 0)
				v = 0;
			if (v > 1)
				v = 1;
			dex = v * (double)len;
			if (dex>len-1)
				dex = len-1;
			*out++ = b[dex];
		}
	}
	else if (channelcount > 1) {
		while (n--) {
			v = *in++;
			if (v < 0)
				v = 0;
			if (v > 1)
				v = 1;
			dex = (long)(v * (double)len) * channelcount;
			if (dex>(len-1)*channelcount)
				dex = (len-1)*channelcount;
			*out++ = b[dex];
		}
	}

	buffer_unlocksamples(buffer);
	return;
zero:
	while (n--)
		*out++ = 0.;
}

// adjust for the appropriate number of inlets and outlets (2 out, no in)
void sgran_dsp64(t_sgran *x, t_object *dsp64, short *count, double samplerate, long maxvectorsize, long flags)
{
	x->w_connected[0] = count[1];
	x->w_connected[1] = count[2];
	object_method(dsp64, gensym("dsp_add64"), x, sgran_perform64, 0, NULL);
}
