# rm_evaluation

---

The module (rm = radiometric measurement) rm_evaluation contains the class RM_Evaluation, which helps with evaluation
of radiometric measurements performed with a radiometric scanning device developed at the Institute of Chemical
Engineering, University Ulm. By creating an instance of the class, it can lead you through the procedure.

Run your python instance in the folder "01Scripts" and put the measurement files in the folder "02Data".

By initializing an instance of the class, the script will scan existing subfolders for exported RM-data. This data
must be stored in a specific folder structure: (Normally created by the export function of the measurement script.)

    folder: 02Data
    |
    |-- folder: date + time + reactor name
    |   |
    |   |-- folder: z-position 1
    |   |   |-- file: *.feather ....................... contains measured spectra in mW nm**-1
    |   |   |-- file: *.log ........................... contains all set parameters for measurement
    |   |   |-- heatmap pictures and *.csv file ....... not used by this script
    |   |
    |   |-- folder: z-position 2
    |   |   |-- ...
    |   |
    |   |-- ...
    |
    |-- ...
    
## Attributes constructed from parameters by \_\_init\_\_:

**saveram** ***(boolean, optional, default = False):*** If True, photon flux spectra are not loaded to use less RAM.

**wl\_ntsr\_borders** ***(tuple of two integers, optional):***
For the calculation of the noise to signal ratio, the spectrum gets separated in two regimes. By default the regime used for defining the noise are the spectrometer pixels 0-150. The signal regime are the pixels 150 to -1
(end of list). This default is defined by the method `import_raw`. Other borders should be defined here as tuple of integers if the automatic import routine is used. The tuple defines the pixel separating the two regimes (predef: 150) and the end of the signal regime.

**predef** ***(list of strings, optional):***
Defines a predefined answer to the import questions of the class when initialized. You give a list of
the answers, it will be red from right to left (`list.pop`).

### Usage

When initializing an instance of the class and `predef == False` you will be asked
if you want to choose a directory to evaluate:

    Choose a Directory to evaluate?
    >?|

Two answers are possible:

- When answering with `'y'` , all export folders in the same directory which have the correct structure and
content will be listed. The program then will ask for the number of the folder
which should be imported.
- If `'n'` is given as answer, all importable .feather files in all subdirectories are listed. The
files can be picked individually by their number.
After choosing the desired data, the program will ask whether it should import these files:

		import raw?
		>?|

If the answer is `'n'` , you can start the import procedure manually with custom parameters by calling
the function `import_raw` . Otherwise the program will also ask if colours should be processed:

    process colors?
    >?|

and for a value of the `ntsr_border` :

    keep ntsr border of 0.550? Or type new
    >?|

After that, all files will be red and processed. Note that this can afford some RAM. The answers to the import questions can be also given in advance by the predef parameter with the class initialization. In case you want to import the second folder ( answer2 = 2 ) in the directory ( answer1 = `'y'` )
and import( answer3 = `'y'` ) without processing the colours ( answer4 = `'n'` ), while keeping the ntsr_border unchanged ( answer5 = `'y'` ), the list `['y', 'n', 'y', 2, 'y']`) must be passed as parameter `predef`. The program reads this list from right to left.

        
## Output Attributes:

**alldfs** ***(numpy array,
        shape: (number of z positions, number of x-pixels, number of y-pixels, number of spectrometer pixels),
        dtype='float32'):***
    Array containing the measured spectra in mW nm<sup>-1</sup>.

**all\_photon\_dfs** ***(numpy array,
        shape: (number of z positions, number of x-pixels, number of y-pixels, number of spectrometer pixels),
        dtype='float32'):***
    Array containing the measured spectra in photons nm<sup>-1</sup> s<sup>-1</sup>.


**allwaves** ***(numpy array,
          shape: (number of z positions, number of spectrometer pixels),
          dtype='float32'):***
    Array containing the respective wavelengths in nm for the spectrometer pixels.

**all\_vars** ***(numpy array,
          shape: (number of z positions),
          dtype='float32'):***
    Array containing a dict with the respective measurement variables:
    
	Name': Name of the Ssan
    
	'z_pos': z-position
    
	'Size': Size of the scan im pixels.
    
	int_time': Integration time.

**all\_integrals** ***(numpy array,
               shape: (number of z positions, number of x-pixels, number of y-pixels),
               dtype='float32'):***
    Containing the received power per pixel in mW.

**all\_photon\_integrals** ***(numpy array,
                      shape: (number of z positions, number of x-pixels, number of y-pixels),
                      dtype='float32'):***
    Containing the received number of photons per pixel.

**all\_rec\_power** ***(numpy array,
              shape: (number of z positions),
              dtype='float32'):***
    Containing the received power per scan in mW.

**all\_rec\_photons** ***(numpy array,
                 shape: (number of z positions),
                 dtype='float32'):***
    Containing the received number of photons per scan.

**all\_nev\_abs** ***(numpy array,
             shape: (number of z positions, number of x-pixels, number of y-pixels, number of spectrometer
             pixels),
             dtype='float32'):***
    Array containing the measured spectra of not absorbable power in mW nm<sup>-1</sup>.

**all\_nev\_abs\_photons** (***numpy array,
                     shape: (number of z positions, number of x-pixels, number of y-pixels, number of
                     spectrometer pixels),
                     dtype='float32'):***
    Array containing the measured spectra of not absorbable photons.

**all\_ratios** ***(numpy array,
            shape: (number of z positions, number of x-pixels, number of y-pixels),
            dtype='float32'):***
    Array containing the peak ratio for every measurement pixel.

**all\_colors** ***(numpy array,
            shape: (number of z positions, number of x-pixels, number of y-pixels, red values, green values,
            blue values),
            dtype='float32'):***
    Array containing the RGB values.

**all\_integrals\_nev\_abs** ***(numpy array,
                       shape: (number of z positions, number of x-pixels, number of y-pixels),
                       dtype='float32'):***
    Containing the not absorbable received power per pixel in mW.

**all\_integrals\_remaining** ***(numpy array,
                         shape: (number of z positions, number of x-pixels, number of y-pixels),
                         dtype='float32'):***
    Containing the not absorbed received power per pixel in mW.

**all\_photon\_integrals\_remaining** ***(numpy array,
                                shape: (number of z positions, number of x-pixels, number of y-pixels),
                                dtype='float32'):***
    Containing the not absorbed received photons per pixel.

## Postprocess and Plot Functions

The following functions are included in the module for postprocesing and plotting. Please refer to the detailed docstrings insde the module's [code](./01Skripts/rm_evaluation.py).


- normalize 
- cos\_corr\_point 
- cos\_corr\_stick 
- show\_histogram 
- SI\_plot 
- import\_raw  
- mak\_peak\_difference\_from\_cache  
- make\_peak\_difference  
- plot 
- plot\_single  
- cutoff\_at\_percent 
- determine\_angle 
- plot\_slice\_set  
- plot\_histograms 
- plot\_slice  
- out\_all 



