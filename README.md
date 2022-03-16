# MUS-177-Final by Angus Yick

This is the final assignment for my music programming class. It is a PD external that shows some techniques such as: lookup table oscillator, wavefolding, allpass filter, amplitude clippers, and tremolo modulator. Main differences from the previous assignment are amplitude clippers and tremolo as well as triangle, ramp, and sawtooth waves. Also, some GUI updates were made overall. Functionality of the external is outlined in the help patch. 

![image](https://user-images.githubusercontent.com/74380180/158534989-2e7a96a4-23cc-47fb-a1d8-c2edcd489bcc.png)

Add the following to the makefile:

In pd_nt:
```
synthesizer~.dll
```

At the bottom:
```
synthesizer~.dll: synthesizer~.c; 
	  cl $(PDNTCFLAGS) $(PDNTINCLUDE) /c $*.c
	  link /dll /export:synthesizer_tilde_setup $*.obj $(PDNTLIB)
```
