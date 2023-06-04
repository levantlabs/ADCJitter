#Created by Manar El-Chammas
#GUI for jitter analysis
#Created Aug. 2022
#Updated May 2023


#Modules
import sys
from PyQt6.QtWidgets import QDialog, QApplication, QPushButton, QRadioButton, QHBoxLayout, QVBoxLayout, QGridLayout, QLineEdit, QLabel
from PyQt6.QtGui import QIntValidator
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from functools import partial
import numpy as np
import pathlib
import time #To save timestamp
import os


#GUI interface
class Window(QDialog):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.debug = False #For debug purposes

        #Custom toolbar definition
        NavigationToolbar.toolitems = (
            ('Home', 'Reset to original view', 'home', 'home'),
            ('Back' , 'Back to previous view', 'back', 'back'),
            ('Forward', 'Forward to next view', 'forward', 'forward'),
            ('Zoom' ,'Zoom to rectangle', 'zoom_to_rect', 'zoom')
        )

        plt.rcParams.update({'font.size': 4})

        self.version = '1.0'

        #Pre-defined texts depending on mode
        self.na_text = 'Not applicable'
        self.clear_text = ''

        #Track number of amplitude & frequency points for custom waveform
        self.ampFreq_numPts = 0

        #Create a figure for plotting
        self.figure = plt.figure()

        #This is the Canvas Widget that displays the `figure`
        self.canvas = FigureCanvas(self.figure)

        #This is the Navigation widget
        self.toolbar = NavigationToolbar(self.canvas, self)

        #Store widget state
        self.app_state = False #If false, then window did not expand to include ampVals

        #Show maximized on startup
        self.showMaximized()


        ################################
        #Define options for signal types
        ################################
        self.choice = 0 #Will be overridden later
        self.choice_button_Label = QLabel(self)
        self.choice_button_Label.setText('Signal Type')

        #Signal options: Sine wave, bandpass, and custom
        self.choice_button1 = QRadioButton('Sine Wave')
        self.choice_button1.setChecked(True)
        self.choice_button1.toggled.connect(self.choiceClicked)

        self.choice_button2 = QRadioButton('Bandpass Signal')
        self.choice_button2.toggled.connect(self.choiceClicked)

        self.choice_button3 = QRadioButton('Custom Signal')
        self.choice_button3.toggled.connect(self.choiceClicked)


        #####################################
        #Additional buttons for functionality
        #####################################
        #Key button: analyzes setup and generates plot
        self.button = QPushButton('Analyze')
        self.button.clicked.connect(self.plot)

        #Saves plot for easy sharing
        self.save = QPushButton('Save')
        self.save.clicked.connect(self.take_screenshot)

        self.title_Label = QLabel(self)
        self.title_Label.setText('First-order Jitter Analysis GUI.  Version {}'.format(self.version))


        ######################################
        #Data entry boxes for analysis
        ######################################

        

        #Minimum frequency
        self.min_frequency_Label = QLabel(self)
        self.min_frequency_Label.setText('Minimum Frequency [MHz]')
        self.min_frequency = QLineEdit(self)
        self.freq_validator_min = QIntValidator(0, 99999, self)
        self.min_frequency.setValidator(self.freq_validator_min)
        #Start with read only
        self.min_frequency.setReadOnly(True)
        self.min_frequency.setText(self.na_text)
        self.min_frequency.textChanged.connect(self.entryValidation_minF)

        #Upper frequency limit (between 0 and high number)
        self.max_Frequency_Label = QLabel(self)
        self.max_Frequency_Label.setText('Maximum Frequency [MHz]')
        self.max_Frequency = QLineEdit(self)
        self.freq_validator_max = QIntValidator(1, 99999, self)
        self.max_Frequency.setValidator(self.freq_validator_max)
        self.max_Frequency.textChanged.connect(self.entryValidation_maxF)

        #Target amplitude
        self.amplitude_Label = QLabel(self)
        self.amplitude_Label.setText('Input Signal Amplitude [dBFS]')
        self.amplitude = QLineEdit(self)
        #self.amplitude.setInputMask("+99999")
        self.ampValidator = QIntValidator(bottom=-400, top=0 , parent=self)
        self.amplitude.setValidator(self.ampValidator)
        self.amplitude.textChanged.connect(self.entryValidation_amp)

        #Set number of frequency amplitude pairs (used in custom signal)
        self.freqResponse_Label = QLabel(self)
        self.freqResponse_Label.setText('Number of Frequency/Amplitude pairs')
        #If the value here changes, call function addFreqResponse
        self.freqResponse = QLineEdit(self)
        response_validator = QIntValidator(2, 9, self)
        self.freqResponse.setValidator(response_validator)
        self.freqResponse.setReadOnly(True)
        self.freqResponse.setText(self.na_text)

        self.freqResponse.textChanged.connect(self.addFreqResponse_2)

         #Amplitude and frequency pairs for custom signal
        self.ampValues = []
        self.freqValues = []
        self.freqValues_label = 'Frequency [MHz]'
        self.ampValues_label = 'Attenuation [dB]'

        self.baseSNR_Label = QLabel(self)
        self.baseSNR_Label.setText('ADC Low Frequency SNR [dBFS]')
        self.baseSNR = QLineEdit(self)

        self.jitter_Label = QLabel(self)
        self.jitter_Label.setText('Target Jitter [psrms]')
        self.jitter = QLineEdit(self)

        self.result_Label = QLabel(self)
        self.result_Label.setText('ADC SNR with above properties [dBFS]')
        self.result = QLineEdit(self)
        self.result.setReadOnly(True)

        self.result2_Label = QLabel(self)
        self.result2_Label.setText('ADC SNR with above properties [dBFS], approximate')
        self.result2 = QLineEdit(self)
        self.result2.setReadOnly(True)

       
        
        
        #Add status text based on run
        self.run_status = QLabel(self)
        self.run_status.setText('Status: Inactive...')
        self.status_text = ''

        # Create the layout
        self.layout1 = QVBoxLayout()
        layout1a = QHBoxLayout()
        
        layout1a.addWidget(self.toolbar)
        layout1a.addWidget(self.save)
        self.layout1.addLayout(layout1a)
        self.layout1.addWidget(self.title_Label)
        self.layout1.addWidget(self.canvas)
        
        layout1b = QHBoxLayout()
        layout2 = QHBoxLayout()
        layout3 = QHBoxLayout()
        layout4 = QHBoxLayout()
        layout5 = QHBoxLayout()
        layout6 = QHBoxLayout()
        layout7 = QHBoxLayout()
        layout_status = QHBoxLayout()
        self.layout8 = QGridLayout()

        layout1b.addWidget(self.choice_button_Label)
        layout1b.addWidget(self.choice_button1)
        layout1b.addWidget(self.choice_button2)
        layout1b.addWidget(self.choice_button3)
        self.layout1.addLayout(layout1b)
        
        #Add input box
        layout2.addWidget(self.min_frequency_Label)
        layout2.addWidget(self.min_frequency)
        layout2.addWidget(self.max_Frequency_Label)
        layout2.addWidget(self.max_Frequency)
        layout3.addWidget(self.amplitude_Label)
        layout3.addWidget(self.amplitude)
        layout3.addWidget(self.freqResponse_Label)
        layout3.addWidget(self.freqResponse)
        

        layout4.addWidget(self.baseSNR_Label)
        layout4.addWidget(self.baseSNR)

        layout4.addWidget(self.jitter_Label)
        layout4.addWidget(self.jitter)

        #layout6.addWidget(self.bw_Label)
        #layout6.addWidget(self.bw)

        layout7.addWidget(self.result_Label)
        layout7.addWidget(self.result)
        layout7.addWidget(self.result2_Label)
        layout7.addWidget(self.result2)

        self.layout1.addLayout(layout2)
        self.layout1.addLayout(layout3)
        #added 5/24/2023
        self.layout1.addLayout(self.layout8)
        self.layout1.addLayout(layout4)
        #self.layout1.addLayout(layout5)
        #layout1.addLayout(layout6)
        
        self.layout1.addWidget(self.button)
        self.layout1.addLayout(layout7)

        #layout7.addWidget(self.attenuation_Label)
        #layout7.addWidget(self.attenuation)

        layout_status.addWidget(self.run_status)

        self.layout1.addLayout(layout_status)


        
        self.setLayout(self.layout1)


    #Entry validation
    def entryValidation_amp(self):
        text = self.amplitude.text()
        self.entryValidation('amplitude', text)
    def entryValidation_maxF(self):
        text = self.max_Frequency.text()
        self.entryValidation('maxFreq', text)
    def entryValidation_minF(self):
        text = self.min_frequency.text()
        self.entryValidation('minFreq', text)
    def entryValidation(self, text, value):
        self.validate(item = text, value=value)
    def validate(self, item, value):
        pos = 0 #This is for the cursor
        self.run_status.setText('Status: Inactive...')

        if value == '':
            return #Nothing to analyze
        if item == 'amplitude':
            result, val, posC = self.ampValidator.validate(value, pos)
            if result == QIntValidator.State.Intermediate  and value != '-':
                if float(val) < self.ampValidator.bottom():
                    self.amplitude.setText(str(self.ampValidator.bottom()))
                elif float(val) > self.ampValidator.top():
                    self.amplitude.setText(str(self.ampValidator.top()))
                    self.status_text = 'Maximum amplitude is 0dBFS.  Reseting to maximum value'
                    self.run_status.setText(self.status_text)
        elif item == 'maxFreq':
            result, val, posC = self.freq_validator_max.validate(value, pos)
            if result == QIntValidator.State.Intermediate  and value != '-':
                if float(val) < self.freq_validator_max.bottom():
                    self.max_Frequency.setText(str(self.freq_validator_max.bottom()))
                    self.status_text = 'Minimum frequency is 1MHz.  Reseting to minimum value'
                    self.run_status.setText(self.status_text)
                elif float(val) > self.freq_validator.top():
                    self.max_Frequency.setText(str(self.freq_validator_max.top()))
        elif item == 'minFreq':
            result, val, posC = self.freq_validator_min.validate(value, pos)
            if result == QIntValidator.State.Intermediate  and value != '-':
                if float(val) < self.freq_validator_min.bottom():
                    self.min_frequency.setText(str(self.freq_validator_min.bottom()))
                    self.status_text = 'Minimum frequency is 1MHz.  Reseting to minimum value'
                    self.run_status.setText(self.status_text)
                elif float(val) > self.freq_validator_min.top():
                    self.min_frequency.setText(str(self.freq_validator_min.top()))
                    
            
                                                 
        
    def choiceClicked(self):
        #If sine wave, gray out min frequency

        if self.choice_button1.isChecked():
            self.min_frequency.setReadOnly(True)
            self.min_frequency.setText(self.na_text)

            self.max_Frequency.setReadOnly(False)
            self.max_Frequency.setText(self.clear_text)

            self.freqResponse.setReadOnly(True)
            self.freqResponse.setText(self.na_text)
            
        #If bandpass, allow text in min and max frequency    
        elif self.choice_button2.isChecked():
            self.min_frequency.setReadOnly(False)
            self.min_frequency.setText(self.clear_text)

            self.max_Frequency.setReadOnly(False)
            self.max_Frequency.setText(self.clear_text)

            self.freqResponse.setReadOnly(True)
            self.freqResponse.setText(self.na_text)
            
        #If custom, gray out min and max frequency
        else:
            self.min_frequency.setReadOnly(True)
            self.min_frequency.setText(self.na_text)

            self.max_Frequency.setReadOnly(True)
            self.max_Frequency.setText(self.na_text)

            self.freqResponse.setReadOnly(False)
            self.freqResponse.setText(self.clear_text)   
        

    def plot(self):
        ''' Plot the defined waveform and results of jitter analysis'''
  
        #Read all the variables based on which waveform is selected
        #Reset status_text, and updated as needed
        
        self.status_text  = ''
        #Sine Wave Analysis
        if self.choice_button1.isChecked():
            
            Fmax = self.max_Frequency.text()
            Fmin = Fmax #self.min_frequency.text()
            self.choice = 1
            self.status_text = 'Sine Wave Analysis. '
        #Bandpass Analysis
        elif self.choice_button2.isChecked():
            Fmax = self.max_Frequency.text()
            Fmin = self.min_frequency.text()
            self.choice = 2
            self.status_text = 'Bandpass signal analysis. '
        #Custom Waveform Analysis
        else:
            self.choice = 3
            Fmax = 0
            Fmin = 0
            self.status_text = 'Custom waveform analysis. '
            if len(self.freqResponse.text()) < 1:
                self.ampFreq_numPts = 0
            else:
                self.ampFreq_numPts = int(self.freqResponse.text())
            

        if self.debug:
            print('\n\n********Choice = {}\n\n'.format(self.choice))

        #Clear figure and create axes and titles
        self.figure.clear()
        ax = self.figure.add_subplot(121)
        ay = self.figure.add_subplot(122)

        ax.title.set_text('Signal Power Spectral Density')
        ay.title.set_text('SNR Curves')

        yL = 'Signal PSD [Normalized]'
        if self.choice == 3:
            yL = 'Signal PSD\n[dB, Normalized]'
           
        ax.set(xlabel='Frequency [MHz]',
               ylabel=yL)
        ay.set(xlabel='Maximum Frequecy [MHz]',
               ylabel='Total SNR [dBFS]')

        
        #Let's start plotting      
        #First, define parameters
        if self.choice == 1: #This is a sine wave.  Only use FMAX
            if len(Fmax) < 1:
                #Set to 100 MHz if no value there
                sigFreq = 100.0
                self.status_text = self.status_text + 'Frequency not set.  Setting to 100MHz. '
            else:
                sigFreq = float(Fmax)
                self.status_text = self.status_text + 'Frequency set to {}MHz. '.format(Fmax)

            bw = 0 #No bandwidth since sine wave
            minF = sigFreq
            maxF = sigFreq

        elif self.choice == 2: #Bandpass signal
            #Make sure Fmax is bigger than Fmin
            if len(Fmin) < 1: #This is empty
                Fmin = 0.0
                if len(Fmax) < 1: #This is also empty
                    Fmax = 100.0
                    self.status_text = self.status_text + 'Frequency not set.  Setting to 0MHz and 100MHz. '
                else:
                    self.status_text = self.status_text + 'Minimum frequency not set.  Setting to 0MHz and {}MHz. '.format(float(Fmax))
            else:
                if len(Fmax) < 1:
                    Fmax = 2*float(Fmin)
                    self.status_text = self.status_text + 'Maximum frequency not set.  Setting to {}MHz and {}MHz. '.format(float(Fmin), Fmax) 
                elif float(Fmax) < float(Fmin): #Fmax is set to less than Fmin
                    self.min_frequency.setText(str(0)) #Force it to 0
                    self.status_text = self.status_text + 'Minimum frequency larger than maximum.  Setting to 0MHz and {}MHz. '.format(float(Fmax))
                    #Should output a status message
                    Fmin = 0.0
                          
            sigFreq = (float(Fmax) + float(Fmin))/2
            bw = float(Fmax) - float(Fmin)
            maxF = sigFreq + bw/2
            minF = sigFreq - bw/2

        elif self.choice == 3: #Custom waveform
            #print(self.ampFreq_numPts)
            if self.ampFreq_numPts == 0: #nothing here, can't analyze
                self.status_text = self.status_text + 'No points entered. Please correct.'
                self.run_status.setText(self.status_text)
                return
            ampVals = np.zeros(self.ampFreq_numPts)
            freqVals = np.zeros(self.ampFreq_numPts)

            #Check if points are valid
            for i in range(self.ampFreq_numPts):
                x1 = self.ampValues[i+1].text()
                x2 = self.freqValues[i+1].text()
                if self.debug:
                    print('Current pair is Freq = {}, Amp = {}'.format(float(x2), float(x1)))

                if len(x1) < 1 or len(x2) < 1:
                    self.status_text = self.status_text + 'Amplitude and frequency points not all defined.  Please correct.'
                    self.run_status.setText(self.status_text)
                    return

                #Now make sure frequency points are in order
                if i < self.ampFreq_numPts - 1:
                    if self.debug:
                        print('Checking frequency = {} vs {}'.format(float(x2), float(self.freqValues[i+2].text())))
                    if float(x2) >= float(self.freqValues[i+2].text()):
                        self.status_text = self.status_text + 'Frequency points are not in ascending order.  Please correct.'
                        self.run_status.setText(self.status_text)
                        return
                       
            for i in range(self.ampFreq_numPts):
                print(self.ampValues[i+1].text())
                #Converting from dB here.  Don't need to convert later
                ampVals[i] = 10**(float(self.ampValues[i+1].text())/20) #+1 because the first one is just a label
                freqVals[i] = float(self.freqValues[i+1].text())

            bw = freqVals[-1] - freqVals[0] #Get full bandwidth
            sigFreq = (freqVals[0] + freqVals[-1]) / 2 #Get center frequency

            maxF = sigFreq + bw/2
            minF = sigFreq - bw/2

        #Get other important parameters
        amp = self.amplitude.text()
        if len(amp) < 1:#Nothing here
            inputAmp = 0 #dBFS
        else:
            inputAmp = float(amp)

        if len(self.baseSNR.text()) < 1: #Nothing here
            adcsnr = 100
            self.status_text = self.status_text + ' Base SNR set to 100dBFS.'
        else:
            adcsnr = float(self.baseSNR.text())

        if len(self.jitter.text()) < 1: #Nothing here
            jitter = 1e-12
            self.status_text = self.status_text + ' Jitter set to 1ps.'
        else:
            jitter = float(self.jitter.text())*1e-12



        #PSD plot is normalized to dBFS (i.e. dBFS doesn't change the y-axis
        #Sine Wave.  Plot a delta function
        if self.choice == 1:
            y = np.asarray([0, 1])
            x = np.asarray([sigFreq, sigFreq])
            #ax.axis([0, (sigFreq)*1.25, 0, 1.25])

        #Bandpass signal
        elif self.choice == 2: #self.ampFreq_numPts == 0: #Assume Rect
            #Plot rect function
            x = np.asarray([minF, minF, maxF, maxF])
            y = np.asarray([0, 1/bw, 1/bw, 0])

        else: #Freq points have been defined
            ampVal_dB = 20*np.log10(ampVals)
            freqPlot = freqVals + 0.0 #This will get plotted in semilogX.  Can't have a 0
            if freqPlot[0] == 0:
                freqPlot[0] = 0.1
            x = np.insert(freqPlot, 0, freqPlot[0])
            x = np.append(x, freqPlot[-1])
            y = np.insert(20*np.log10(ampVals), 0, 0)
            y = np.insert(ampVal_dB, 0, min(ampVal_dB)-2)
            y = np.append(y, min(ampVal_dB)-2)

        if self.choice == 3: #Plot with log axis for custom waveform
            ax.semilogx(x, y, 'k')
        else:
            ax.plot(x, y, 'k')
        if self.choice == 1:
            ax.plot(sigFreq, 1, 'k^') #Add an arrow for fun
        ax.axis([0, maxF*1.25, min(0, min(y)), 1.25*max(y)])


        #All required parameters for jitter analysis are available.
        #Can now calculate impact of jitter

        #Normalize signal power to 1 for simplicity
        Ps_unit = 0.5 #Amplitude is 1 --> power is 1^2/2 = 0.5
        #Convert power to the amplitude dBFS
        Ps = Ps_unit * 10**(inputAmp/10)
        if self.debug:
            print('Signal Power = {}'.format(Ps))

        baseNoise = Ps_unit / 10**(adcsnr / 10)
        if self.debug:
            print('Base Noise Power = {}'.format(baseNoise))

        sineAmp = np.sqrt(Ps*2)
        if self.debug:
            print('Actual sine amplitude = {}'.format(sineAmp))

        #Calculate SNR across several different jitter values
        rg = 4
        lgt = []
        maxSNR = 0
        minSNR = 1000000
        jitter_array = np.zeros(rg)
        snr_array = np.zeros(rg)
        for j in range(rg):
            jit = jitter / 2.0 * float(1+j)
            jitter_array[j] = jit
                                       
            numPts = 20
            freqRange = np.geomspace(1, sigFreq*1.5, numPts)
            #freqRange = sigFreq

            #For sine wave analysis
            if self.choice == 1:
                #Plot across range of frequencies
                jitterNoise = 2*np.pi**2*(freqRange*1e6)**2*sineAmp**2*jit**2
                #jitterNoise = np.zeros(len(freqRange))
                #for f in range(len(freqRange)):   
                #    jitterNoise[0], snrtemp = self.jitter_approximate(sigFreq, bw, Ps, Ps_unit, jit, baseNoise)
                #jitterActual = 2*np.pi**2 * (sigFreq*1e6)**2*sineAmp**2*jitter**2
                jitterActual = (2*np.pi* (sigFreq*1e6) * sineAmp)**2 / 2 * jitter**2
            #Bandpass signal
            elif self.choice == 2: 
                #For now, assuming rect from minF to maxF
                #Geom space can't include 0
                freqRange = np.geomspace(max(minF, 1), maxF*1.5, numPts)
                #sigFreq = maxF
                if minF == 0:
                    jitterNoise = 2*np.pi**2*(freqRange*1e6)**2/3*sineAmp**2*jit**2
                    jitterActual = 2*np.pi**2 * (maxF*1e6)**2/3*sineAmp**2*jitter**2
                    #jitterActual = (2*np.pi * (maxF * 1e6)
                else:
                    
                    #freqRange = np.geomspace(bw, maxF*2, numPts)
                    #print(freqRange)
                    jitterNoise = np.zeros(len(freqRange))
                    for f in range(len(freqRange)):
                        centerF = (freqRange[f] - minF)*0.5 + minF
                        bw_temp = freqRange[f] - minF
                        jitterNoise[f], snrtemp = self.jitter_approximate(centerF, bw_temp, Ps, Ps_unit, jit, baseNoise)
                        jitterActual = 0

            elif self.choice == 3: 
                #freqRange = np.geomspace(bw, maxF*2, numPts)
                freqFactor = np.geomspace(0.1, 2, numPts) #Scale bandwidth
                freqRange = freqFactor * maxF
                jitterNoise = np.zeros(len(freqRange))
                
                if self.debug:
                    print('Amplitude in voltage = {}'.format(ampVals))
                
                for f in range(len(freqFactor)):
                    
                    jitterNoise[f], snrtemp = self.jitter_approximate(0, 0, Ps, Ps_unit, jit, baseNoise, ampVal = ampVals, freqVal = freqVals * freqFactor[f])
                    jitterActual = 0

            #Calculate in dBFS
            snr_ratio = Ps_unit / (baseNoise + jitterNoise) #Simulated
            if self.debug:
                if self.choice == 1:
                    print('PS_UNIT = {}, jitterActual = {}'.format(Ps_unit, jitterActual))
            snr_ratio_actual = Ps_unit / (baseNoise + jitterActual) #Theoretical
            realSNR = 10*np.log10(snr_ratio)
            #if len(realSNR) > 1:
            maxSNR_temp = max(realSNR)
            minSNR_temp = min(realSNR)
            #else:
            #    maxSNR_temp = realSNR
            #    minSNR_temp = minSNR
                
            if maxSNR_temp > maxSNR:
                maxSNR = maxSNR_temp
            if minSNR_temp < minSNR:
                minSNR = minSNR_temp
                
            snrActual = 10*np.log10(snr_ratio_actual)

            if self.choice == 1: #bw == 0:
                self.result.setText('{0:.2f}'.format(snrActual))
            elif self.choice == 2:
                if minF == 0: #Theoretical estimate
                    self.result.setText('{0:.2f}'.format(snrActual))
                else:
                    self.result.setText('N/A')
            elif self.choice == 3:
                self.result.setText('N/A')
        
            ay.semilogx(freqRange, realSNR)
            
            lgt.append('{0:.2f} psrms'.format(jit/1e-12))

            #Add line at actual target frequency
            #Calculate approximation based on properties
            if self.choice == 1:
                ampVals = []
                freqVals = []
            elif self.choice == 2: #self.ampFreq_numPts == 0:
                minF_estimate = sigFreq - bw/2
                maxF_estimate = sigFreq + bw/2
                Npts = (maxF_estimate - minF_estimate) / 1 + 1 #In steps of 1MHz
                freqPts = np.arange(Npts) / (Npts-1) * (maxF_estimate - minF_estimate) + minF_estimate
                freqVals = freqPts
              
                ampVals = np.zeros(len(freqPts))
                for i in range(int(Npts)):
                    ampVals[i] = 1 #Rect function

            

        if self.debug:
            if self.choice == 1:
                print('Parameters for SNR Estimate.  sigFreq = {}, bw = {}, Ps = {}, Ps_unit = {}, jitter = {}, baseNoise = {}'.format(sigFreq, bw, Ps, Ps_unit, jitter, baseNoise))
                
        jitter_estimate, snr_estimate_dB = self.jitter_approximate(sigFreq, bw, Ps, Ps_unit, jitter, baseNoise, ampVal = ampVals, freqVal = freqVals)
        self.result2.setText('{0:0.2f}'.format(snr_estimate_dB))
                
           

        ay.legend(lgt)
        ay.semilogx([maxF, maxF], [minSNR, maxSNR], 'k-.')

        plt.tight_layout()
        #ax.tighten_layout()

        # refresh canvas
        self.canvas.draw()
        

        self.run_status.setText(self.status_text)

        
    def jitter_approximate(self, sigFreq, bw, Ps, Ps_unit, jitter, baseNoise, ampVal = [], freqVal = []):
        #Add line at actual target frequency
        #Calculate approximation based on properties
        #amp and freq are for transfer functions.  If empty, do the stuff below
        
  
        #Custom waveform
        if len(freqVal) > 0 and len(ampVal) > 0: #Then, work on this data and create a full spectrum
            Npts = 100
            Fstep = (max(freqVal) - min(freqVal)) / Npts
            freqPts = np.zeros(Npts)
            ampPts = np.zeros(Npts)
            if freqVal[0] == 0: #0 complicates the math without value
                freqVal[0] = 0.01
            #ampVal = ampVal ** 2 #Set to power here, and then interpolate in power
            #Get a finer definition of the amplitude and the frequency
            for i in range(Npts):
                currFreq = freqVal[0] + i*Fstep #Start from the first frequency
                #Find the frequency point to the left of currFreq
                ind = np.searchsorted(freqVal, currFreq, side='left')
                #print(ind)
                #Interpolate the amplitude
                #Linear interpolation
                #currAmp = (ampVal[ind] - ampVal[ind-1]) / (freqVal[ind] - freqVal[ind-1]) * (currFreq - freqVal[ind-1]) + ampVal[ind-1]
                lowFreq = freqVal[ind-1]
                #if freqVal[ind-1] == 0:
                #    lowFreq = 0.001
                #else:
                #    lowFreq = freqVal[ind-1]
                #if currFreq == 0:
                #    currFreq = 0.001
                if self.debug:
                    print('*** Low Freq = {}, ampVal[ind] = {}, ampVal[ind-1] = {}, freqVal[ind] = {}, currFreq = {}'.format(lowFreq, ampVal[ind], ampVal[ind-1], freqVal[ind], currFreq))
                currAmp_dB = (20*np.log10(ampVal[ind]) - 20*np.log10(ampVal[ind-1])) * np.log10(currFreq / lowFreq) / np.log10(freqVal[ind] / lowFreq) + 20*np.log10(ampVal[ind-1])
                #Working through the math
                #10** --> currFreq / freqVal
                currAmp = 10**(currAmp_dB / 20)
                freqPts[i] = currFreq
                ampPts[i] = currAmp

            amp = ampPts
        
        else: #if len(freqVal) == 0:
            minF_estimate = sigFreq - bw/2
            maxF_estimate = sigFreq + bw/2
            Fstep = 1# 0.1
            Npts = int((maxF_estimate - minF_estimate) / Fstep + 1) #In steps of 0.1MHz
            freqPts = np.arange(Npts) / (Npts) * (maxF_estimate - minF_estimate) + minF_estimate

        #elif len(ampVal) == 0:
            amp = np.zeros(len(freqPts))
            for i in range(int(Npts)):
                amp[i] = 1 #Rect function


        #Normalize total transfer function to 1 (np.sum(amp)) and then scale to power
        
        amp = np.sqrt(1 * amp / np.sum(amp) * Ps/Ps_unit)#Normalized
        
        if self.debug:
            print('Amplitude = {}'.format(amp))
            print('Signal power ideal = {}, estimated = {}'.format(Ps, np.sum(amp**2/2)))

        #Now, go through and integrate jitter noise on this estimate
        jitter_estimate = 0
        for i in range(int(Npts)):
            #Calculate jitter for every slice of the bandwidth
            jitter_estimate +=  2*np.pi**2 * (freqPts[i]*1e6)**2*amp[i]**2 *jitter**2
        if self.debug:
            print('Jitter_estimate = {}'.format(jitter_estimate))

        #In dBFS, so use full-scale power
        snr_estimate = Ps_unit / (baseNoise + jitter_estimate)
        snr_estimate_dB = 10*np.log10(snr_estimate)

        return jitter_estimate, snr_estimate_dB
        
    def take_screenshot(self):
        pixmap = self.grab()
        if not pixmap.isNull():
            if self.debug:
                print('Pixmap is not null')
            primary_path = str(pathlib.Path().resolve())
            sigType = ''
            if self.choice == 1:
                sigType = 'sine_'
            elif self.choice == 2:
                sigType = 'bandpass_'
            elif self.choice == 3:
                sigType = 'custom_'
            image_dir = primary_path + '/images'
            image_path = image_dir + '/plot_' + sigType + str(time.time()) + '.png'

            if not os.path.exists(image_dir):
                os.makedirs(image_dir)
            
            
            result = pixmap.save(image_path, "PNG")
            if result:
                print('Image saved successfully...')
            else:
                print('Image not saved...')
            self.run_status.setText('Image saved to {}'.format(image_path))
        else:
            self.run_status('Error: Image not saved.')
        #QApplication.quit()

          

    def addFreqResponse_2(self): #Second attempt
        if len(self.freqResponse.text()) < 1 or self.freqResponse.text() == self.na_text: #Don't do anything
               return

        numPts = int(self.freqResponse.text())
        #print('Number of points in frequency response = {}'.format(numPts))
        #currLen = len(self.ampValues) #The first one is just the label
        #print('Length of ampValues = {}'.format(currLen))
    
        if numPts > 0 and len(self.ampValues) < 1: #add label
            self.ampValues.append(QLabel(self))
            self.freqValues.append(QLabel(self))
            self.ampValues[0].setText(self.ampValues_label)
            self.freqValues[0].setText(self.freqValues_label)
            #self.ampValues.append(self.ampValues_label)
            #self.freqValues.append(self.freqValues_label)
            self.layout8.addWidget(self.ampValues[0], 0, 0)
            self.layout8.addWidget(self.freqValues[0], 1, 0)
            

        #Now, add in all the others
        currLen = len(self.ampValues)-1 #The first one is just the label

        #Two cases.  Case 1: numPts > currLen.  Case 2: numPts < currLen
        if numPts > currLen: #Remove header
            #Add in points
            for i in range(numPts):
                if 1+i > currLen:
                    self.ampValues.append(QLineEdit(self))
                    self.freqValues.append(QLineEdit(self))
                    self.layout8.addWidget(self.ampValues[-1], 0, 1+i)
                    self.layout8.addWidget(self.freqValues[-1], 1, 1+i)
        elif numPts < currLen:
            #Remove points
            for i in range(currLen):
                if 1+i > numPts:
                     self.layout8.removeWidget(self.ampValues[-1])
                     self.layout8.removeWidget(self.freqValues[-1])
                     #Below from stack overflow.  Purge widgets properly
                     self.ampValues[-1].deleteLater()
                     self.freqValues[-1].deleteLater()
                     self.ampValues.pop()
                     self.freqValues.pop()
                   

        if numPts == 0 and len(self.ampValues) > 0: #Remove the widgets
            print('Removing heading...')
            print('Length of ampValues = {}'.format(len(self.ampValues)))
            self.layout8.removeWidget(self.ampValues[0])
            self.layout8.removeWidget(self.freqValues[0])
            self.ampValues[0].deleteLater()
            self.freqValues[0].deleteLater()
            self.ampValues.pop()
            self.freqValues.pop()
            


        return
               
        



    
    
#Main program
if __name__ == '__main__':
    #Define application and show window
    app = QApplication(sys.argv)
    main = Window()
    main.show()
    sys.exit(app.exec())
