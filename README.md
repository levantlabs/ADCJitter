# ADC Jitter Analysis GUI

Clock jitter limits ADC performance.  The impact of jitter depends on the signal statistics.  This GUI provides a simple interface to analyze this impact.

## Why This Matters

As signal frequencies increase, ADC performance becomes more sensitive to clock jitter.  The SNR degradation is a function of both the input signal statistics and the jitter variance.  

Jitter analysis is typically done with sine wave input signals.  However, wideband signals can relax the jitter specification, allowing for lower power designs.  

## How To Use It

Thie GUI allows for three different signal types:
- Sine wave
- Bandpass power spectral density (PSD)
- A custom power spectral density

The signal type is selected with the radio button at the top of the gui.

Additional parameters are defined, in order to calculate the impact.  This includes:
- Minimum frequency: this is the minimum frequency in the PSD.  This is only applicable in the bandpass signal type
- Maximum frequency:  this is the maximum frequency in the PSD.  This is only applicable with a sine wave signal and the bandpass signal type
- Number of frequency/amplitude pairs:  this is only applicable for the custom waveform.  When a non-zero number, two additional rows are added to the gui
  - Attenuation & Frequency:  This are the attenuation and frequency pairs that define the custom waveform
- ADC low frequency SNR:  this is the base SNR for the ADC without jitter
- Target jitter: this is the target jitter for the analysis

Many of these do not need values.  If the box is left empty, the tool will fill in the required values.  For example, the target jitter is set to 1psrms, and the maximum frequency is set to 100MHz.

When all the parameters are entered, click on the **Analyze** button.  This will calculate the SNR is possible, and will also approximate the SNR.  Two plots are generated.  The first is the signal PSD.  The second is the impact of frequency on the iitter curve.

Finally, the screen can be saved by clicking on the save button.

## Setting The Tool Up

A `requirements.txt` file is included, and can be used by pip to set up a virtual environment with all the required libraries.

The tool is opened by running ``python main.py``.

## Final Note

This tool is under construction.  The calulated jitter for the custom waveform is still under evaluation.


